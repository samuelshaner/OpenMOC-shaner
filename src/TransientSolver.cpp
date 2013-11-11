#include "TransientSolver.h"

/**
 * Transient constructor
 * @param geom pointer to the geometry
 * @param cmfd pointer to the cmfdTransient
 * @param solver pointer to the solver
 */
TransientSolver::TransientSolver(Geometry* geom, CmfdTransient* cmfd, Solver* solver){

  log_printf(NORMAL, "Creating transient object...");
  
  _cmfd = cmfd;
  _geom = geom;
  _solver = solver;
  _mesh = _cmfd->getMesh();
  _materials = _mesh->getMaterials();
  _solve_method = _mesh->getSolveType();
  _num_groups = _geom->getNumEnergyGroups();
  _timer = new Timer();
 
  /* compute the volume of the core */
  computeVolCore();

  _mesh->createNewFlux(REFERENCE);

}

/**
 * Transient Destructor clears all memory
 */
TransientSolver::~TransientSolver() {
}


void TransientSolver::solveInitialState(){

    log_printf(NORMAL, "Solving time dependent problem...");
    
    /* initialize the time stepper and transient material properties */
    initializeTimeStepper();
    initializeTransientMaterials();
    
    /* compute initial shape function */
    _k_eff_0 = static_cast<Cmfd*>(_cmfd)->computeKeff();
    
    /* set the initial eigenvalue in the cmfd solver */
    log_printf(NORMAL, "k_eff_0: %f", _k_eff_0);
    _cmfd->setKeff0(_k_eff_0);
    
    _mesh->copyFlux(CURRENT, PREVIOUS);

    /* compute the power normalization factor */
    _power_factor = _power_init / computePower();
    _cmfd->vecScale(_mesh->getFluxes(PREVIOUS), _power_factor);
    _cmfd->vecScale(_mesh->getFluxes(CURRENT), _power_factor);
        
    /* initialize the precursor concentration and PKE matrices */
    initializePrecursorConc();
    copyPrecConc(CURRENT, PREVIOUS);
        
    _cmfd->setInitialState(false);
    
    /* save initial global parameters */
    _time.push_back(0.0);
    _temp.push_back(computeCoreTemp());
    _power.push_back(computePower());
    log_printf(NORMAL, "TIME: %f, POWER: %.10f, TEMP: %f", _time.back(), _power.back(), _temp.back());  
}


void TransientSolver::solveOuterStep(){
  
    /* set time integration parameters */
    double tolerance_out = 1.e-6;
    int length_conv = _power.size();
    double residual_out = 1.0;
    double keff;
    int iter = 0;
    
    if (_ts->getTime(CURRENT) == _start_time)
	_timer->startTimer();
    
    /* compute the flux at the forward time step */
    _ts->takeStep();
    syncMaterials(CURRENT);
    updatePrecursorConc();    

    while (residual_out > tolerance_out){
	_mesh->copyFlux(CURRENT, REFERENCE);
	_cmfd->computeKeff();
	updateTemperatures();
	updatePrecursorConc();
	syncMaterials(CURRENT);	
	residual_out = computeResidual();
	iter++;
    }
        
    /* compute and save the power, temp, and time */
    _time.push_back(_ts->getTime(CURRENT));
    _temp.push_back(computeCoreTemp());
    _power.push_back(computePower());
    log_printf(NORMAL, "TIME: %f, POWER: %.10f, TEMP: %f, ITER: %i", _time.back(), _power.back(), _temp.back(), iter);  
        
    copyFieldVariables(CURRENT, PREVIOUS);

    if (fabs(_ts->getTime(CURRENT) - _end_time) < 1e-8){
	std::string msg_string;
	log_printf(TITLE, "TIMING REPORT");
	_timer->stopTimer();
	_timer->recordSplit("Total time to solve transient problem");
	
	double tot_time = _timer->getSplit("Total time to solve transient problem");
	msg_string = "Total time to solve transient problem";
	msg_string.resize(53, '.');
	log_printf(RESULT, "%s%4.4E sec", msg_string.c_str(), tot_time);
    }    
}


double TransientSolver::computeResidual(){
    
    double res = 0.0;
    
    double* new_flux = _mesh->getFluxes(CURRENT);
    double* fwd_flux = _mesh->getFluxes(REFERENCE);
    
    for (int i = 0; i < _mesh->getNumCells(); i++){
	for (int e = 0; e < _num_groups; e++)
	    res += pow(_materials[i]->getNuSigmaF()[e] * (fwd_flux[e] - new_flux[e]) * _mesh->getVolumes()[i]/ (_materials[i]->getNuSigmaF()[e] * fwd_flux[e] * _mesh->getVolumes()[i] + 1e-10), 2.0);
    }
    
    res = pow(res, 0.5);
    res = res / (_mesh->getNumCells() * _num_groups);
    
    return res;
}

void TransientSolver::copyFieldVariables(materialState state_from, materialState state_to){

    copyTemperature(state_from, state_to);
    _mesh->copyFlux(state_from, state_to);
    copyPrecConc(state_from, state_to);
    
}


void TransientSolver::copyTemperature(materialState state_from, materialState state_to){
    
    /* loop over mesh cells */
    for (int i = 0; i < _mesh->getNumCells(); i++)
	_materials[i]->copyTemperature(state_from, state_to); 
}


void TransientSolver::trimVectors(int len_conv){

    int len_curr = _power.size();
    for (int i = len_conv; i < len_curr; i++){
	_time.pop_back();
	_power.pop_back();
	_temp.pop_back();
    }
}


/* change xs due to temperature changes and move control blade */
void TransientSolver::updateTemperatures(){

    log_printf(INFO, "Performing temperature feedback...");
    
    /* initialize variables */
    double power, temp;
    double* flux = _mesh->getFluxes(CURRENT);
    
    /* loop over fsrs */
    for (int i = 0; i < _mesh->getNumCells(); i++){
	
	power = 0.0;
	
	/* compute power in fsr */
	for (int e = 0; e < _num_groups; e++){
	    power += _alpha / _nu * _materials[i]->getNuSigmaF()[e] * flux[i*_num_groups + e];
	}
		
	power = power;
	temp = _materials[i]->getTemperature(PREVIOUS) + _dt_cmfd * power;
	_materials[i]->setTemperature(CURRENT, temp);
    }
}



/* compute core temperature */
double TransientSolver::computeCoreTemp(){

    log_printf(INFO, "Computing core temperature...");
    
    double geom_temp = 0.0;
    
    /* loop over mesh cells */
    for (int i = 0; i < _mesh->getNumCells(); i++){    
	if (_materials[i]->isFissionable()){
	    geom_temp += _materials[i]->getTemperature(CURRENT) * _mesh->getVolumes()[i];
	}
    }
    
    geom_temp = geom_temp / _vol_core;
    
    return geom_temp;
}


void TransientSolver::copyPrecConc(materialState state_from, materialState state_to){
    
    /* loop over mesh cells */
    for (int i = 0; i < _mesh->getNumCells(); i++){    
	if (_materials[i]->getType() == FUNCTIONAL){
	    static_cast<FunctionalMaterial*>(_materials[i])->copyPrecConc(state_from, state_to);
	}
    } 
}



/* initialize the precursors for each cell */
void TransientSolver::initializePrecursorConc(){
  
    log_printf(INFO, "Initializing precursors...");
    
    double power = 0.0;
    
    double* flux = _mesh->getFluxes(CURRENT);
    
    for (int i = 0; i < _mesh->getNumCells(); i++){
	
	for (int dg = 0; dg < _num_delay_groups; dg++){
	    power = 0.0;
	    
	    for (int e = 0; e < _num_groups; e++){
		power += _materials[i]->getNuSigmaF()[e] * flux[i*_num_groups + e] / _k_eff_0;
	    }
	    
	    if (power > 0.0){
		static_cast<FunctionalMaterial*>(_materials[i])->setPrecConc(REFERENCE,     _cmfd->getBeta()[dg] / _cmfd->getLambda()[dg] * power, dg);
		
		log_printf(DEBUG, "cell: %i, dg: %i, initial prec: %f", i, dg, _cmfd->getBeta()[dg] / _cmfd->getLambda()[dg] * power);
	    }
	}
    }

    copyPrecConc(REFERENCE, PREVIOUS);
    copyPrecConc(REFERENCE, PREVIOUS_CONV);
    copyPrecConc(REFERENCE, CURRENT);
    copyPrecConc(REFERENCE, FORWARD);
}


void TransientSolver::updatePrecursorConc(){

    log_printf(INFO, "Updating precursor concentrations...");
    
    /* initialize variables */
    double power, old_conc, new_conc;

    double* flux = _mesh->getFluxes(CURRENT);

    /* loop over fsrs */
    for (int i = 0; i < _mesh->getNumCells(); i++){
	
	for (int dg = 0; dg < _num_delay_groups; dg++){
	    
	    power = 0.0;
	    
	    /* compute power in fsr */
	    for (int e = 0; e < _num_groups; e++)
		power += _materials[i]->getNuSigmaF()[e]*flux[i*_num_groups+e] / _k_eff_0;	    
	    
	    if (_materials[i]->isFissionable()){
		old_conc = static_cast<FunctionalMaterial*>(_materials[i])->getPrecConc(PREVIOUS, dg);
		new_conc = old_conc / (1.0 + _cmfd->getLambda()[dg]*_dt_cmfd) 
		    + _cmfd->getBeta()[dg]*_dt_cmfd*power/(1.0 + _cmfd->getLambda()[dg]*_dt_cmfd);
		
		log_printf(DEBUG, "cell: %i, old conc: %f, new conc: %f", i, old_conc, new_conc);
		
		static_cast<FunctionalMaterial*>(_materials[i])->setPrecConc(CURRENT, new_conc, dg);
	    }
	}
    }
}


/* compute the power / cc */
double TransientSolver::computePower(materialState state){

    log_printf(INFO,"Computing power...");
    
    double power = 0.0;
    double* flux = _mesh->getFluxes(state);
    double* volumes = _mesh->getVolumes();
    
    for (int i = 0; i < _mesh->getNumCells(); i++){
	for (int e = 0; e < _num_groups; e++)
	    power += _kappa / _nu * _materials[i]->getNuSigmaF()[e]*flux[i*_num_groups + e] * volumes[i];
    }
    
    log_printf(DEBUG, "core power: %f, core power density: %f", power, power / _vol_core);
    
    power = power / _vol_core;
    
    return power;
}


void TransientSolver::setStartTime(double time){
    _start_time = time;
}


void TransientSolver::setEndTime(double time){
    _end_time = time;
}


void TransientSolver::initializeTimeStepper(){
    _ts = new TimeStepper(_start_time, _end_time, _dt_moc, _dt_cmfd);
    _cmfd->initialize(_ts);
}


void TransientSolver::initializeTransientMaterials(){

    for (int i = 0; i < _mesh->getNumCells(); i++){
	if (_materials[i]->getType() == FUNCTIONAL){
	    static_cast<FunctionalMaterial*>(_materials[i])->setTimeStepper(_ts);
	    static_cast<FunctionalMaterial*>(_materials[i])->initializeTransientProps(_num_delay_groups);
	}
  }
}


void TransientSolver::setDtMOC(double dt){
    _dt_moc = dt;
}


void TransientSolver::setDtCMFD(double dt){
    _dt_cmfd = dt;
    _cmfd->setDtCMFD(dt);
}


void TransientSolver::computeVolCore(){

    _vol_core = 0.0;
    
    for (int i = 0; i < _mesh->getNumCells(); i++){
	
	if (_materials[i]->isFissionable()){
	    log_printf(DEBUG, "cell: %i, vol: %f", i, _mesh->getVolumes()[i]);
	    _vol_core += _mesh->getVolumes()[i];
	}
    }
}


void TransientSolver::setKappa(double kappa){
    _kappa = kappa;
}


void TransientSolver::setNu(double nu){
    _nu = nu;
}


void TransientSolver::setAlpha(double alpha){
    _alpha = alpha;
}


void TransientSolver::setNumDelayGroups(int num_groups){
    _num_delay_groups = num_groups;
    _cmfd->setNumDelayGroups(num_groups);
}


void TransientSolver::setTransientMethod(const char* trans_type){

    if (strcmp("ADIABATIC", trans_type) == 0)
	_transient_method = ADIABATIC;
    else
	log_printf(ERROR, "Could not recognize transient method; "
		   " the options are ADIABATIC");
    
    _cmfd->setTransientType(_transient_method);
}


void TransientSolver::setPowerInit(double power){
    _power_init = power;
}


void TransientSolver::syncMaterials(materialState state){

    for (int i = 0; i < _mesh->getNumCells(); i++){
	if (_materials[i]->getType() == FUNCTIONAL)
	    static_cast<FunctionalMaterial*>(_materials[i])->sync(state); 
    }
}
