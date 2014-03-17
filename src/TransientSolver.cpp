#include "TransientSolver.h"

/**
 * Transient constructor
 * @param geom pointer to the geometry
 * @param cmfd pointer to the cmfdTransient
 * @param solver pointer to the solver
 */
TransientSolver::TransientSolver(Geometry* geom, Tcmfd* tcmfd, Cmfd* cmfd, Solver* solver){

    log_printf(NORMAL, "Creating transient object...");
    
    _tcmfd = tcmfd;
    _cmfd = cmfd;
    _geom = geom;
    _solver = solver;
    _mesh = _geom->getMesh();
    _geom_mesh = _geom->getGeomMesh();
    _solve_method = _mesh->getSolveType();
    _ng = _geom->getNumEnergyGroups();
    _timer = new Timer();
    _temp_peak_value = 0.0;
    _power_peak_value = 0.0;
    _amp_solve_counter = 0;
    _prolongation = true;
    
    /* compute the volume of the core */
    computeVolCore();

    _mesh->createNewFlux(PREVIOUS_CONV);    
    _mesh->createNewFlux(PREVIOUS);
    _mesh->createNewFlux(CURRENT);
    _mesh->createNewFlux(FORWARD);
    _mesh->createNewFlux(FORWARD_PREV);

    /* create new flux arrays */
    if (_solve_method == MOC){
	_geom_mesh->createNewFlux(FORWARD);
	_geom_mesh->createNewFlux(CURRENT);
	_geom_mesh->createNewFlux(PREVIOUS_CONV);
	_geom_mesh->createNewFlux(PREVIOUS);
	_geom_mesh->eraseFlux(SHAPE);
	_geom_mesh->setFlux(SHAPE, _mesh->getFSRFluxes());
    }
}

/**
 * Transient Destructor clears all memory
 */
TransientSolver::~TransientSolver() {
}


void TransientSolver::solveInitialState(){

    log_printf(NORMAL, "Solving time dependent problem...");

    /* start timer for global solve */
    _timer->startTimer();

    /* initialize the time stepper and transient material properties */
    initializeTimeStepper();
    initializeTransientMaterials();
    initializeTransientLogfile();
    syncMaterials(CURRENT);

    /* compute initial shape function */
    log_printf(NORMAL, "Computing initial shape");
    if (_solve_method == DIFFUSION){
	_k_eff_0 = _cmfd->computeKeff();
	_mesh->setKeff0(_k_eff_0);
	_power_factor = _power_init / computePower();
	vecScale(_mesh->getFluxes(CURRENT), _power_factor, _mesh->getNumCells()*_ng);
	_mesh->copyFlux(CURRENT, PREVIOUS);
    }
    else{

	/* MOC solve */
        static_cast<ThreadPrivateSolverTransient*>(_solver)->resetSegmentMaterials();        
	_k_eff_0 = static_cast<ThreadPrivateSolverTransient*>(_solver)->convergeSource(1000);    
	_mesh->setKeff0(_k_eff_0);

	/* compute fine mesh flux shape */
	_mesh->copyFlux(SHAPE, CURRENT);
	_mesh->computeFineShape(_geom_mesh->getFluxes(CURRENT), _mesh->getFluxes(CURRENT));

	/* normalize vectors to initial power level */
	_power_factor = _power_init / computePower();
	vecScale(_mesh->getFluxes(CURRENT), _power_factor, _mesh->getNumCells()*_ng);
	vecScale(_geom_mesh->getFluxes(SHAPE), _power_factor, _geom_mesh->getNumCells()*_ng);

	/* copy flux to FORWARD, PREVIOUS, and PREVIOUS_CONV */
	_mesh->copyFlux(CURRENT, SHAPE);
	_mesh->copyFlux(CURRENT, FORWARD);
	_mesh->copyFlux(CURRENT, PREVIOUS);
	_mesh->copyFlux(CURRENT, PREVIOUS_CONV);

	/* copy NDCs to CURRENT and PREVIOUS_CONV */
	_mesh->copyDs(FORWARD, CURRENT);
	_mesh->copyDs(FORWARD, PREVIOUS_CONV);

	/* copy fine mesh shape to PREVIOUS, PREVIOUS_CONV, and FORWARD */
	_geom_mesh->copyFlux(CURRENT, PREVIOUS);
	_geom_mesh->copyFlux(CURRENT, PREVIOUS_CONV);
	_geom_mesh->copyFlux(CURRENT, FORWARD);
    } 

    /* set the initial eigenvalue in the cmfd solver */
    log_printf(NORMAL, "k_eff_0: %f", _k_eff_0);
    _tcmfd->setKeff0(_k_eff_0);

    /* initialize the precursor concentration and PKE matrices */
    initializePrecursorConc();

    if (_solve_method == MOC){
	mapPrecConc();
	_tcmfd->computeFrequency();
    }

    /* set initial state flag to false */
    _tcmfd->setInitialState(false);
    
    /* save initial global parameters */
    _num_amp_solves.push_back(0);
    _time_tcmfd.push_back(0.0);
    _time_moc.push_back(0.0);
    _temp_core.push_back(computeCoreTemp());
    _power_core.push_back(computePower());
    _power_peak.push_back(_power_peak_value);
    _temp_peak.push_back(_temp_peak_value);
    _reactivity.push_back(static_cast<ThreadPrivateSolverTransient*>(_solver)->getReactivity());
    _moc_iters.push_back(static_cast<ThreadPrivateSolverTransient*>(_solver)->getNumIters());
    logStep();
    log_printf(NORMAL, "TIME: %f, POWER: %.10f, TEMP: %f", 
	       _time_tcmfd.back(), _power_core.back(), _temp_core.back());  
}


void TransientSolver::solveOuterStep(){
  
    /* set time integration parameters */
    double tolerance_in = 1.e-8;
    double tolerance_out = 1.e-4;
    double residual_out = 1.0;
    double residual_in = 1.0;
    int iter = 0;
    double k_eff;
    int iter_outer = 0;
    double old_power = _power_core.back();
    
    if (_solve_method == MOC){
	
	/* reset outer loop residual */
	residual_out = 1.0;
	int len_conv = _power_core.size();
	updatePrecursorConc(PREVIOUS, FORWARD);

	/* converge the source on the outer time step fine mesh */
	while (residual_out > tolerance_out){

	    if (iter != 0 || !_prolongation){
	        /* compute the shape at the forward time step */
	        sync(FORWARD);
		k_eff = static_cast<ThreadPrivateSolverTransient*>(_solver)->convergeSource(1000);
		_mesh->computeFineShape(_geom_mesh->getFluxes(FORWARD), _mesh->getFluxes(CURRENT));

		/* if ADIABATIC method, renormalize coarse mesh flux */
		if (_transient_method == ADIABATIC){
		    _geom_mesh->copyFlux(FORWARD, CURRENT);
		    vecScale(_geom_mesh->getFluxes(FORWARD), _power_core.back() / computePower(), 
			     _ng*_geom_mesh->getNumCells()); 
		}
	    }
	    
	    /* set FORWARD_PREV coarse flux to test for convergence */
	    _mesh->copyFlux(CURRENT, FORWARD_PREV);

	    /* trim power, temp, and time to account due to implicit solve */
	    trimVectors(len_conv);
	        
	    /* reset CURRENT time */
	    _ts->setTime(CURRENT, _ts->getTime(PREVIOUS_CONV));
	    _ts->setTime(PREVIOUS, _ts->getTime(PREVIOUS_CONV));

	    /* reset precursor conc, temp, and flux to PREVIOUS_CONV */
	    copyFieldVariables(PREVIOUS_CONV, CURRENT);
	    copyFieldVariables(PREVIOUS_CONV, PREVIOUS);
	    _mesh->copyFlux(PREVIOUS_CONV, CURRENT);
	    _mesh->copyFlux(PREVIOUS_CONV, PREVIOUS);

	    /* inner loop over TCMFD solves */
	    while (_ts->getTime(CURRENT) < _ts->getTime(PREVIOUS_CONV) + _ts->getDtMOC() - 1e-8){
		
		/* reset iterator and residual */
		iter = 0;
		residual_in = 1.0;
		
		/* increment CURRENT time */
		_ts->takeStep();
		
		/* sync cross sections and interpolate flux/Ds */
		sync(CURRENT);

		/* update precursor concentration */
		updatePrecursorConc();

		/* take outer step */
		while (residual_in > tolerance_in){
		    _mesh->copyFlux(CURRENT, FORWARD);
		    _timer->startTimer();
		    _tcmfd->solveTCMFD();
		    _timer->stopTimer();
		    _timer->recordSplit("TCMFD solve");
		    updateTemperatures();
		    updatePrecursorConc();
		    sync(CURRENT);
		    residual_in = computeResidual(DIFFUSION);
		    iter++;
		}
		
		/* copy the coarse mesh flux from CURRENT to PREVIOUS */
		_mesh->copyFlux(CURRENT, PREVIOUS);
		copyFieldVariables(CURRENT, PREVIOUS);
	
		/* compute and save the power, temp, and time */
		_num_amp_solves.push_back(iter);
		_time_tcmfd.push_back(_ts->getTime(CURRENT));
		_temp_core.push_back(computeCoreTemp());
		_power_core.push_back(computePower());
		_power_peak.push_back(_power_peak_value);
		_temp_peak.push_back(_temp_peak_value);
		log_printf(NORMAL, "TIME: %f, POWER: %.10f, TEMP: %f, ITER: %i", _time_tcmfd.back(), _power_core.back(), _temp_core.back(), iter);  		
	    }

	    /* compute frequencies */
	    _tcmfd->computeFrequency();

	    /* compute and print residual */
	    residual_out = computeResidual(MOC);
	    log_printf(NORMAL, "tolerance: %1.3E, resid out: %1.3E", tolerance_out, residual_out);

	    iter_outer++;

	    /* output to logfile */
	    _reactivity.push_back(static_cast<ThreadPrivateSolverTransient*>(_solver)->getReactivity());
	    _moc_iters.push_back(static_cast<ThreadPrivateSolverTransient*>(_solver)->getNumIters());
	    logStep();

	    /* update fine mesh flux a forward step */
	    _geom_mesh->reconstructFineFlux(_geom_mesh->getFluxes(FORWARD), _mesh->getFluxes(CURRENT));

	}

	/* copy converged field variables */
	copyFieldVariables(CURRENT, PREVIOUS_CONV);
	copyFieldVariables(CURRENT, FORWARD);
	_mesh->copyFlux(CURRENT, PREVIOUS_CONV);
	_mesh->copyDs(FORWARD, PREVIOUS_CONV);

	/* set PREVIOUS_CONV time */
	_ts->convergedMOCStep();

	log_printf(NORMAL, "Iterations to converge shape function: %i", iter_outer);
    }
    else{
	while (_ts->getTime(CURRENT) < _ts->getTime(PREVIOUS_CONV) + _ts->getDtMOC() - 1e-8){
	    iter = 0;
	    _ts->takeStep();
	    sync(CURRENT);
	    updatePrecursorConc();
	    residual_in = 1.0;

	    /* take outer step */
	    while (residual_in > tolerance_in){
		_mesh->copyFlux(CURRENT, FORWARD);
		_timer->startTimer();
		_tcmfd->solveTCMFD();
		_timer->stopTimer();
		_timer->recordSplit("TCMFD solve");
		updateTemperatures();
		updatePrecursorConc();
		sync(CURRENT);
		residual_in = computeResidual(DIFFUSION);
		iter++;
	    }
	  
	    copyFieldVariables(CURRENT, PREVIOUS);
  
	    /* compute and save the power, temp, and time */
	    _time_tcmfd.push_back(_ts->getTime(CURRENT));
	    _temp_core.push_back(computeCoreTemp());
	    _power_core.push_back(computePower());
	    log_printf(NORMAL, "TIME: %f, POWER: %.10f, TEMP: %f, ITER: %i", _time_tcmfd.back(), _power_core.back(), _temp_core.back(), iter);  	    
	}
	
	_ts->convergedMOCStep();
    }
    
    if (fabs(_ts->getTime(CURRENT) - _end_time) < 1e-8){
	std::string msg_string;
	log_printf(TITLE, "TIMING REPORT");
	_timer->stopTimer();
	_timer->recordSplit("Total time to solve transient problem");
        _timer->printSplits();
	
	double tot_time = _timer->getSplit("Total time to solve transient problem");
	msg_string = "Total time to solve transient problem";
	msg_string.resize(53, '.');
	log_printf(RESULT, "%s%4.4E sec", msg_string.c_str(), tot_time);
    }    
}


double TransientSolver::computeResidual(solveType solve_type){

    _timer->startTimer();
    
    double res = 0.0;
    double source = 0.0;
    double power_new = 0.0;
    double power_old = 0.0;
    
    if (solve_type == MOC){

      
	std::vector<int>::iterator iter;
	double* fine_shape = _geom_mesh->getFluxes(FORWARD);
	double* coarse_flux_new = _mesh->getFluxes(CURRENT);
	double* coarse_flux_old = _mesh->getFluxes(FORWARD_PREV);
	double* volumes = _mesh->getVolumes();
	double* velocity = _mesh->getVelocity();
	Material** materials = _geom_mesh->getMaterials();
	
	for (int i = 0; i < _mesh->getNumCells(); i++){
	    
	    for (iter = _mesh->getCellFSRs()->at(i).begin(); iter != _mesh->getCellFSRs()->at(i).end(); ++iter){

		for (int e = 0; e < _ng; e++){
		    if (materials[*iter]->isFissionable()){
		        power_new = materials[*iter]->getNuSigmaF()[e] * coarse_flux_new[i*_ng+e] * volumes[i] * fine_shape[*iter*_ng+e] / velocity[e]; 
			power_old = materials[*iter]->getNuSigmaF()[e] * coarse_flux_old[i*_ng+e] * volumes[i] * fine_shape[*iter*_ng+e] / velocity[e]; 
				  
			res += pow((power_new - power_old) / power_new, 2.0);
		    }
		}
	    }
	}

	res = res / (_geom_mesh->getNumCells() * _ng);
	res = pow(res, 0.5);
      
	/*
      for (int i = 0; i < _mesh->getNumCells(); i++){
	for (int e = 0; e < _ng; e++){
	  if (_mesh->getFluxes(CURRENT)[i*_ng+e] != 0.0)
	    res += pow((_mesh->getFluxes(CURRENT)[i*_ng+e] - _mesh->getFluxes(FORWARD_PREV)[i*_ng+e]) / _mesh->getFluxes(CURRENT)[i*_ng+e], 2);
	}
      }
    
      res = pow(res, 0.5);
      res = res / (_mesh->getNumCells() * _ng);
	*/

    }
    else{
	double* new_flux = _mesh->getFluxes(CURRENT);
	double* old_flux = _mesh->getFluxes(FORWARD);
	Material** materials = _mesh->getMaterials();
	double* volumes = _mesh->getVolumes();

	
	for (int i = 0; i < _mesh->getNumCells(); i++){
	    for (int e = 0; e < _ng; e++){
		if (materials[i]->isFissionable())
		    res += pow((materials[i]->getNuSigmaF()[e] * 
				(new_flux[i*_ng+e] - old_flux[i*_ng+e]) * volumes[i]) 
			       / (materials[i]->getNuSigmaF()[e] * new_flux[i*_ng+e] * volumes[i]), 2.0);
	    }
	}

	res = pow(res, 0.5);	
	res = res / (_mesh->getNumCells() * _ng);

    }
    
    _timer->stopTimer();
    _timer->recordSplit("Computing residual");
    
    return res;
}

void TransientSolver::copyFieldVariables(materialState state_from, materialState state_to){
  
    _timer->startTimer();

    copyTemperature(state_from, state_to);
    copyPrecConc(state_from, state_to);
    
    if (_solve_method == MOC)
	_geom_mesh->copyFlux(state_from, state_to);
    else
	_mesh->copyFlux(state_from, state_to);

    _timer->stopTimer();
    _timer->recordSplit("Copying field variables");
}


void TransientSolver::copyTemperature(materialState state_from, materialState state_to){

    if (_solve_method == MOC){
	for (int i = 0; i < _geom_mesh->getNumCells(); i++)
	    _geom_mesh->getMaterials()[i]->copyTemperature(state_from, state_to); 
    }
    else{
	for (int i = 0; i < _mesh->getNumCells(); i++)
	    _mesh->getMaterials()[i]->copyTemperature(state_from, state_to); 	
    }    
}


void TransientSolver::trimVectors(int len_conv){

    int len_curr = _power_core.size();
    for (int i = len_conv; i < len_curr; i++){
	_time_tcmfd.pop_back();
	_power_core.pop_back();
	_temp_core.pop_back();
	_temp_peak.pop_back();
	_power_peak.pop_back();
	_num_amp_solves.pop_back();
    }
}


/* change xs due to temperature changes and move control blade */
void TransientSolver::updateTemperatures(){

    log_printf(INFO, "Performing temperature feedback...");

    _timer->startTimer();
    
    /* initialize variables */
    double power_old, power_new, power, temp;

    if (_solve_method == MOC){
	std::vector<int>::iterator iter;
	double* fine_shape_new = _geom_mesh->getFluxes(CURRENT);
	double* coarse_flux_new = _mesh->getFluxes(CURRENT);
	double* fine_shape_old = _geom_mesh->getFluxes(PREVIOUS);
	double* coarse_flux_old = _mesh->getFluxes(PREVIOUS);
	double* volumes = _mesh->getVolumes();
	double* velocity = _mesh->getVelocity();
	Material** materials = _geom_mesh->getMaterials();
	
	/* loop over coarse mesh cells */
        #pragma omp parallel for private(power_old, power_new, power, temp, iter)
	for (int i = 0; i < _mesh->getNumCells(); i++){
	    
	    /* loop over fsrs within cell */
	    for (iter = _mesh->getCellFSRs()->at(i).begin(); iter != _mesh->getCellFSRs()->at(i).end(); ++iter){

		power_new = 0.0;
		power_old = 0.0;

		/* compute power in fsr */
		for (int e = 0; e < _ng; e++){
		    power_new += _alpha * materials[*iter]->getSigmaF()[e] 
			* fine_shape_new[*iter*_ng+e] * coarse_flux_new[i*_ng+e] 
		      * volumes[i] / velocity[e];

		    power_old += _alpha * materials[*iter]->getSigmaF()[e] 
			* fine_shape_old[*iter*_ng+e] * coarse_flux_old[i*_ng+e] 
		      * volumes[i] / velocity[e];
		    
		}
		
		power = (power_old + power_new) / 2.0;
		temp = materials[*iter]->getTemperature(PREVIOUS) + _dt_cmfd * power;
		materials[*iter]->setTemperature(CURRENT, temp);		
	    }
	}
    }
    else{
	double* flux_new = _mesh->getFluxes(CURRENT);
	double* flux_old = _mesh->getFluxes(PREVIOUS);
	Material** materials = _mesh->getMaterials();
	
	/* loop over coarse mesh cells */
        #pragma omp parallel for private(power_old, power_new, power, temp)
	for (int i = 0; i < _mesh->getNumCells(); i++){
	    
	    power_old = 0.0;
	    power_new = 0.0;
	    
	    /* compute power in fsr */
	    for (int e = 0; e < _ng; e++){
		power_new += _alpha * materials[i]->getSigmaF()[e] * flux_new[i*_ng+e];
		power_old += _alpha * materials[i]->getSigmaF()[e] * flux_old[i*_ng+e];
	    }
	    
	    power = (power_old + power_new) / 2.0;
	    temp = materials[i]->getTemperature(PREVIOUS) + _dt_cmfd * power;
	    materials[i]->setTemperature(CURRENT, temp);		
	}	
    }

    _timer->stopTimer();
    _timer->recordSplit("Updating temperature");
}



/* compute core temperature */
double TransientSolver::computeCoreTemp(){

    log_printf(INFO, "Computing core temperature...");

    _timer->startTimer();
    
    double geom_temp = 0.0;
    _temp_peak_value = 0.0;

    if (_solve_method == MOC){
	/* loop over mesh cells */
	for (int i = 0; i < _geom_mesh->getNumCells(); i++){    
	    if (_geom_mesh->getMaterials()[i]->isFissionable()){
		geom_temp += _geom_mesh->getMaterials()[i]->getTemperature(CURRENT) * _geom_mesh->getVolumes()[i];
		_temp_peak_value = std::max(_temp_peak_value, _geom_mesh->getMaterials()[i]->getTemperature(CURRENT));
	    }
	}
    }
    else{
	/* loop over mesh cells */
	for (int i = 0; i < _mesh->getNumCells(); i++){    
	    if (_mesh->getMaterials()[i]->isFissionable()){
	        geom_temp += _mesh->getMaterials()[i]->getTemperature(CURRENT) * _mesh->getVolumes()[i];
		_temp_peak_value = std::max(_temp_peak_value, _mesh->getMaterials()[i]->getTemperature(CURRENT));
	    }
	}
    }
        
    geom_temp = geom_temp / _vol_core;
    
    _timer->stopTimer();
    _timer->recordSplit("Computing core temp");

    return geom_temp;
}


void TransientSolver::copyPrecConc(materialState state_from, materialState state_to){

    if (_solve_method == MOC){
	/* loop over mesh cells */
        #pragma omp parallel for
	for (int i = 0; i < _geom_mesh->getNumCells(); i++){    
	    if (_geom_mesh->getMaterials()[i]->getType() == FUNCTIONAL){
		static_cast<FunctionalMaterial*>(_geom_mesh->getMaterials()[i])->copyPrecConc(state_from, state_to);
	    }
	} 	
    }
    else{
	/* loop over mesh cells */
        #pragma omp parallel for
	for (int i = 0; i < _mesh->getNumCells(); i++){    
	    if (_mesh->getMaterials()[i]->getType() == FUNCTIONAL){
		static_cast<FunctionalMaterial*>(_mesh->getMaterials()[i])->copyPrecConc(state_from, state_to);
	    }
	} 	
    }
}



/* initialize the precursors for each cell */
void TransientSolver::initializePrecursorConc(){
  
    log_printf(INFO, "Initializing precursors...");
    
    double power = 0.0;    

    if (_solve_method == MOC){
	
	FP_PRECISION* flux = _geom_mesh->getFluxes(SHAPE);
	Material** materials = _geom_mesh->getMaterials();

	for (int i = 0; i < _geom_mesh->getNumCells(); i++){
	    
	    for (int dg = 0; dg < _ndg; dg++){
		power = 0.0;

		for (int e = 0; e < _ng; e++){
		    power += materials[i]->getNuSigmaF()[e] * flux[i*_ng + e] / _k_eff_0;
		}

		if (power > 0.0){
		    static_cast<FunctionalMaterial*>(materials[i])->setPrecConc(CURRENT, _tcmfd->getBeta()[dg] / _tcmfd->getLambda()[dg] * power, dg);
		  
		    log_printf(DEBUG, "cell: %i, dg: %i, initial prec: %f", i, dg, _tcmfd->getBeta()[dg] / _tcmfd->getLambda()[dg] * power);
		}
	    }
	}

	copyPrecConc(CURRENT, PREVIOUS_CONV);
	copyPrecConc(CURRENT, PREVIOUS);
    }
    else{
	double* flux = _mesh->getFluxes(CURRENT);
	Material** materials = _mesh->getMaterials();

	for (int i = 0; i < _mesh->getNumCells(); i++){
	    
	    for (int dg = 0; dg < _ndg; dg++){
		power = 0.0;
		
		for (int e = 0; e < _ng; e++){
		    power += materials[i]->getNuSigmaF()[e] * flux[i*_ng + e] / _k_eff_0;
		}
		
		if (power > 0.0){
		    static_cast<FunctionalMaterial*>(materials[i])->setPrecConc(CURRENT, _tcmfd->getBeta()[dg] / _tcmfd->getLambda()[dg] * power, dg);
		    
		    log_printf(DEBUG, "cell: %i, dg: %i, initial prec: %f", i, dg, _tcmfd->getBeta()[dg] / _tcmfd->getLambda()[dg] * power);
		}
	    }
	}
	
	copyPrecConc(CURRENT, PREVIOUS);
	copyPrecConc(CURRENT, PREVIOUS_CONV);
    }    
}


void TransientSolver::updatePrecursorConc(materialState state_from, materialState state_to){

    log_printf(INFO, "Updating precursor concentrations...");

    _timer->startTimer();
    
    /* initialize variables */
    double power_new, power_old, old_conc, new_conc, k1, k2, k3;   
    
    if (_solve_method == MOC){
	std::vector<int>::iterator iter;
	double* fine_shape_new = _geom_mesh->getFluxes(state_to);
	double* fine_shape_old = _geom_mesh->getFluxes(state_from);
	double* coarse_flux_new = _mesh->getFluxes(state_to);
	double* coarse_flux_old = _mesh->getFluxes(state_from);
	double* volumes = _mesh->getVolumes();
	double* velocity = _mesh->getVelocity();
	Material** materials = _geom_mesh->getMaterials();
	double* lambda = _tcmfd->getLambda();
	double* beta = _tcmfd->getBeta();
	double dt = _ts->getTime(state_to) - _ts->getTime(state_from);
	
	/* loop over coarse mesh cells */
        #pragma omp parallel for private(power_new, power_old, old_conc, new_conc, k1, k2, k3, iter)
	for (int i = 0; i < _mesh->getNumCells(); i++){

	    /* loop over fsrs within cell */
	    for (iter = _mesh->getCellFSRs()->at(i).begin(); iter != _mesh->getCellFSRs()->at(i).end(); ++iter){
		
		/* loop over delayed groups */
		for (int dg = 0; dg < _ndg; dg++){
		    
		    power_old = 0.0;
		    power_new = 0.0;

		    /* compute ks */
		    k1 = exp(- lambda[dg] * dt);
		    k2 = 1.0 - (1.0 - k1) / (lambda[dg] * dt);
		    k3 = k1  - (1.0 - k1) / (lambda[dg] * dt);

		    /* compute power in fsr */
		    for (int e = 0; e < _ng; e++){
			power_new += materials[*iter]->getNuSigmaF()[e] * fine_shape_new[*iter*_ng+e] 
			    * coarse_flux_new[i*_ng+e] * volumes[i] / velocity[e] / _k_eff_0;
			power_old += materials[*iter]->getNuSigmaF()[e] * fine_shape_old[*iter*_ng+e] 
			    * coarse_flux_old[i*_ng+e] * volumes[i] / velocity[e] / _k_eff_0;
		    }		    

		    if (materials[*iter]->isFissionable()){
			old_conc = static_cast<FunctionalMaterial*>(materials[*iter])->getPrecConc(state_from, dg);
			new_conc = k1 * old_conc + k2 * beta[dg] / lambda[dg] * power_new - k3 * beta[dg] / lambda[dg] * power_old;
			
			log_printf(DEBUG, "fsr: %i, old conc: %.12f, new conc: %.12f", *iter, old_conc, new_conc);
			static_cast<FunctionalMaterial*>(materials[*iter])->setPrecConc(CURRENT, new_conc, dg);
		    }
		}
	    }
	}

	mapPrecConc();
    }
    else{
	double* flux_new = _mesh->getFluxes(CURRENT);
	double* flux_old = _mesh->getFluxes(PREVIOUS);
	Material** materials = _mesh->getMaterials();
	double* lambda = _tcmfd->getLambda();
	double* beta = _tcmfd->getBeta();
	
	/* loop over mesh cells */
	#pragma omp parallel for private(power_new, power_old, old_conc, new_conc, k1, k2, k3)
	for (int i = 0; i < _mesh->getNumCells(); i++){
	    
	    for (int dg = 0; dg < _ndg; dg++){
		
		power_old = 0.0;
		power_new = 0.0;

		/* compute ks */
		k1 = exp(- lambda[dg] * _dt_cmfd);
		k2 = 1.0 - (1.0 - k1) / (lambda[dg] * _dt_cmfd);
		k3 = k1  - (1.0 - k1) / (lambda[dg] * _dt_cmfd);
		
		/* compute power in fsr */
		for (int e = 0; e < _ng; e++){
		    power_new += materials[i]->getNuSigmaF()[e]*flux_new[i*_ng+e] / _k_eff_0;
		    power_old += materials[i]->getNuSigmaF()[e]*flux_old[i*_ng+e] / _k_eff_0;
		
		}

		if (materials[i]->isFissionable()){
		    old_conc = static_cast<FunctionalMaterial*>(materials[i])->getPrecConc(PREVIOUS, dg);
		    new_conc = k1 * old_conc + k2 * beta[dg] / lambda[dg] * power_new - k3 * beta[dg] / lambda[dg] * power_old;
		    
		    log_printf(DEBUG, "cell: %i, old conc: %f, new conc: %f", i, old_conc, new_conc);
		    
		    static_cast<FunctionalMaterial*>(materials[i])->setPrecConc(CURRENT, new_conc, dg);
		}
	    }
	}	
    }

    _timer->stopTimer();
    _timer->recordSplit("Updating precursor concentration");
}


/* map the precursor concentration fron the geometry mesh to the CMFD mesh */
void TransientSolver::mapPrecConc(){

    _timer->startTimer();

    double prec_conc;
    
    /* interator to loop over fsrs in each mesh cell */
    std::vector<int>::iterator iter;
    Material** fine_materials = _geom_mesh->getMaterials();
    double* fine_volumes = _geom_mesh->getVolumes();
    Material** coarse_materials = _mesh->getMaterials();
    double* coarse_volumes = _mesh->getVolumes();

    /* loop over CMFD cells */
    #pragma omp parallel for private(prec_conc, iter)
    for (int i = 0; i < _mesh->getNumCells(); i++){
	
	for (int dg = 0; dg < _ndg; dg++){
	    
	    prec_conc = 0.0;
	    
	    /* loop over FSRs in mesh cell */
	    for (iter = _mesh->getCellFSRs()->at(i).begin(); iter != _mesh->getCellFSRs()->at(i).end(); ++iter){
		if (fine_materials[*iter]->getType() == FUNCTIONAL){
		    prec_conc += static_cast<FunctionalMaterial*>(fine_materials[*iter])->getPrecConc(CURRENT, dg) * fine_volumes[*iter];
		}
	    }	    
	    
	    if (coarse_materials[i]->getType() == FUNCTIONAL){
		static_cast<FunctionalMaterial*>(coarse_materials[i])->setPrecConc(CURRENT, prec_conc / coarse_volumes[i], dg);
	    }
	}
    }    

    _timer->stopTimer();
    _timer->recordSplit("Mapping precursor concentration");
}


/* compute the power / cc */
double TransientSolver::computePower(materialState state){

    log_printf(INFO, "Computing power...");

    _timer->startTimer();
    
    double power = 0.0;
 
    if (_solve_method == MOC){
	std::vector<int>::iterator iter;
	double* fine_shape = _geom_mesh->getFluxes(CURRENT);
	double* coarse_flux = _mesh->getFluxes(CURRENT);
	double* volumes = _mesh->getVolumes();
	double* fine_volumes = _geom_mesh->getVolumes();
	double* velocity = _mesh->getVelocity();
	Material** materials = _geom_mesh->getMaterials();
	double power_cell;
	_power_peak_value = 0.0;
	
	for (int i = 0; i < _mesh->getNumCells(); i++){
	    for (iter = _mesh->getCellFSRs()->at(i).begin(); 
		 iter != _mesh->getCellFSRs()->at(i).end(); ++iter){	 
	        power_cell = 0.0;
	        for (int e = 0; e < _ng; e++){
		    power += _kappa * materials[*iter]->getSigmaF()[e]
			* fine_shape[*iter*_ng + e] * coarse_flux[i*_ng+e] * volumes[i] 
			/ velocity[e] * fine_volumes[*iter];
		    power_cell += _kappa * materials[*iter]->getSigmaF()[e]
			* fine_shape[*iter*_ng + e] * coarse_flux[i*_ng+e] 
			/ velocity[e] * fine_volumes[*iter];
		}

		_power_peak_value = std::max(_power_peak_value, power_cell); 
	    }
	}	    
    }
    else{
	double* flux = _mesh->getFluxes(state);
	double* volumes = _mesh->getVolumes();
	
	for (int i = 0; i < _mesh->getNumCells(); i++){
	    for (int e = 0; e < _ng; e++){
	        power += _kappa * _mesh->getMaterials()[i]->getSigmaF()[e]
		  * flux[i*_ng + e] * volumes[i];
		_power_peak_value = std::max(_power_peak_value, _kappa * flux[i*_ng + e]
					     * _mesh->getMaterials()[i]->getSigmaF()[e]);  
	    }
	}
    }
    
    log_printf(DEBUG, "core power: %f, core power density: %f", power, power / _vol_core);	
    
    power = power / _vol_core;

    _timer->stopTimer();
    _timer->recordSplit("Computing power");
    
    return power;
}


void TransientSolver::setStartTime(double time){
    _start_time = time;
}


void TransientSolver::setEndTime(double time){
    _end_time = time;
}


void TransientSolver::initializeTimeStepper(){
    log_printf(NORMAL, "initializing time stepper");
    _ts = new TimeStepper(_start_time, _end_time, _dt_moc, _dt_cmfd);
    _tcmfd->setTimeStepper(_ts);
    _mesh->setTimeStepper(_ts);
}


void TransientSolver::initializeTransientMaterials(){

    log_printf(NORMAL, "initializing transient materials...");

    if (_solve_method == MOC){
	for (int i = 0; i < _geom_mesh->getNumCells(); i++){
	    /* clone materials */
	    if (_geom_mesh->getMaterials()[i]->getType() == BASE){
		_geom_mesh->getMaterials()[i] = _geom_mesh->getMaterials()[i]->clone();
	    }
	    else{
		_geom_mesh->getMaterials()[i] = static_cast<FunctionalMaterial*>(_geom_mesh->getMaterials()[i])->clone();
		static_cast<FunctionalMaterial*>(_geom_mesh->getMaterials()[i])->setTimeStepper(_ts);
		static_cast<FunctionalMaterial*>(_geom_mesh->getMaterials()[i])->initializeTransientProps(_ndg, false);
	    }
	}
    }

    for (int i = 0; i < _mesh->getNumCells(); i++){
	if (_mesh->getMaterials()[i]->getType() == FUNCTIONAL){
	    static_cast<FunctionalMaterial*>(_mesh->getMaterials()[i])->setTimeStepper(_ts);
	    static_cast<FunctionalMaterial*>(_mesh->getMaterials()[i])->initializeTransientProps(_ndg, true);
	}
    }
    
    log_printf(NORMAL, "done initializing transient materials...");
}


void TransientSolver::initializeTransientLogfile(){

  /* set logfile name */
  _log_file = std::string(getLogfileName());
  _log_file.resize(_log_file.length()-4);
  _log_file.append("-transient.log");

  /* open file */
  std::ofstream logfile;
  logfile.open((std::string(getOutputDirectory()) + "/" + _log_file).c_str(), std::ios::app);
  logfile << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << std::endl;
  logfile << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << std::endl;
  logfile << "- - - - - - - - - TRANSIENT OUTPUT FILE - - - - - - - - - -" << std::endl;
  logfile << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << std::endl;
  logfile << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << std::endl;
  
  /* write trial summary */
  logfile << "Number of fine mesh cells   = " << _geom_mesh->getNumCells() << std::endl;
  logfile << "Number of coarse mesh cells = " << _mesh->getNumCells() << std::endl;

  if (_solve_method == DIFFUSION)
    logfile << "Shape fct method            = DIFFUSION" << std::endl;
  else{
    logfile << "Shape fct method            = MOC" << std::endl;
    if (_transient_method == MAF)
      logfile << "Transient method            = MAF" << std::endl;
    else if (_transient_method == ADIABATIC)
      logfile << "Transient method            = ADIABATIC" << std::endl;
    if (_transient_method == THETA)
      logfile << "Transient method            = THETA" << std::endl;
  }
  
  logfile << "Time step size inner        = " << _dt_cmfd << " s " << std::endl;
  logfile << "Time step size outer        = " << _dt_moc << " s " << std::endl;

  /* write file heading */
  logfile << std::endl;
  logfile << "Time (s)\tCore Power\tPeak Power\tCore Temp (C)\tPeak Temp (C)\t# Amp Solves\t# MOC iters\tReactivity ($)" << std::endl;
  
  logfile.close();
}


void TransientSolver::logStep(){

  /* open file */
  std::ofstream logfile;
  logfile.open((std::string(getOutputDirectory()) + "/" + _log_file).c_str(), std::ios::app);
  logfile.precision(5);
  logfile << std::scientific;
  int num_steps = round(_dt_moc / _dt_cmfd);
  int i = 0;

  if (_num_amp_solves.size() > 1)
    i = _num_amp_solves.size() - num_steps;

  if (_num_amp_solves.size() > 1){
    for (i = _num_amp_solves.size() - num_steps; i < _num_amp_solves.size()-1; i++){
      logfile << _time_tcmfd.at(i) << "\t" << _power_core.at(i) << "\t" << _power_peak.at(i) 
	      << "\t"  << _temp_core.at(i) << "\t" << _temp_peak.at(i) << "\t" << _num_amp_solves.at(i)
	      << "\t\t-\t\t-\t\t" << std::endl;
    }
  }

  /* write final line */
  logfile << _time_tcmfd.at(i) << "\t" << _power_core.at(i) << "\t" << _power_peak.at(i) 
	  << "\t" << _temp_core.at(i) << "\t" << _temp_peak.at(i) << "\t" << _num_amp_solves.at(i) 
	  << "\t\t" << _moc_iters.back() << "\t\t" << _reactivity.back() << std::endl;

  /* close file */
  logfile.close();
}


void TransientSolver::setDtMOC(double dt){
    _dt_moc = dt;
    _mesh->setDtMOC(dt);
    _geom_mesh->setDtMOC(dt);

    if (_ts != NULL)
      _ts->setDtMOC(dt);

}


void TransientSolver::setDtCMFD(double dt){
    _dt_cmfd = dt;
    _tcmfd->setDtCMFD(dt);

    if (_ts != NULL)
      _ts->setDtCMFD(dt);

}


void TransientSolver::computeVolCore(){

    _vol_core = 0.0;

    if (_solve_method == MOC){
	for (int i = 0; i < _geom_mesh->getNumCells(); i++){
	    
	    if (_geom_mesh->getMaterials()[i]->isFissionable()){
		log_printf(DEBUG, "fsr: %i, vol: %f", i, _geom_mesh->getVolumes()[i]);
		_vol_core += _geom_mesh->getVolumes()[i];
	    }
	}
    }
    else{
	for (int i = 0; i < _mesh->getNumCells(); i++){
	    
	    if (_mesh->getMaterials()[i]->isFissionable()){
		log_printf(DEBUG, "cell: %i, vol: %f", i, _mesh->getVolumes()[i]);
		_vol_core += _mesh->getVolumes()[i];
	    }
	}	
    }    
    
    log_printf(INFO, "core volume: %f", _vol_core);
}


void TransientSolver::setKappa(double kappa){
    _kappa = kappa;
}


void TransientSolver::setAlpha(double alpha){
    _alpha = alpha;
}


void TransientSolver::setNumDelayGroups(int num_groups){
    _ndg = num_groups;
    _tcmfd->setNumDelayGroups(num_groups);
    _mesh->setNumDelayGroups(num_groups);
    _geom_mesh->setNumDelayGroups(num_groups);
}


void TransientSolver::setTransientMethod(const char* trans_type){
 
    if (strcmp("ADIABATIC", trans_type) == 0){
	_transient_method = ADIABATIC; 
	log_printf(NORMAL, "Transient method set to ADIABATIC");
    }
    else if (strcmp("THETA", trans_type) == 0){
	_transient_method = THETA;
	log_printf(NORMAL, "Transient method set to THETA");
    }
    else if (strcmp("MAF", trans_type) == 0){
	_transient_method = MAF;
	log_printf(NORMAL, "Transient method set to MAF");
    }
    else
	log_printf(ERROR, "Could not recognize transient method; "
		   " the options are ADIABATIC, THETA, and MAF");

    _tcmfd->setTransientType(_transient_method);
}


void TransientSolver::setPowerInit(double power){
    _power_init = power;
}


void TransientSolver::syncMaterials(materialState state){

   if (_solve_method == MOC){

        /* sync the fine mesh material cross sections */
        #pragma omp parallel for
	for (int i = 0; i < _geom_mesh->getNumCells(); i++){
	    if (_geom_mesh->getMaterials()[i]->getType() == FUNCTIONAL){
		static_cast<FunctionalMaterial*>(_geom_mesh->getMaterials()[i])->sync(state); 
	    }
	}
    }
    else{
        #pragma omp parallel for
	for (int i = 0; i < _mesh->getNumCells(); i++){
	    if (_mesh->getMaterials()[i]->getType() == FUNCTIONAL){
		static_cast<FunctionalMaterial*>(_mesh->getMaterials()[i])->sync(state); 
	    }
	}	
    }


}


void TransientSolver::sync(materialState state){

    _timer->startTimer();

    if (_solve_method == MOC){

        /* sync the fine mesh material cross sections */
        syncMaterials(state);	

	/* interpolate the fine mesh shape */
	_geom_mesh->interpolateFlux(_ts->getImprovedRatio(state));

	/* reconstruct fine mesh flux */
	_mesh->reconstructFineFlux(_geom_mesh->getFluxes(CURRENT), _mesh->getFluxes(state)); 

	/* compute coarse mesh XS */
	_mesh->computeXS(_geom_mesh, state);

	/* copy coarse flux to CURRENT */
	_mesh->copyFlux(state, CURRENT);

	/* compute course mesh d_hats */
	_mesh->computeDs(1.0, CURRENT);

	/* interpolate NDCs */
	_mesh->interpolateDs(_ts->getImprovedRatio(state));
    }
    else{
        syncMaterials(state);
    }

    _timer->stopTimer();
    _timer->recordSplit("Syncing");
}


double TransientSolver::getPower(){
    return _power_core.back();
}


double TransientSolver::getTemp(){
    return _temp_core.back();
}


double TransientSolver::getTime(){
    return _time_tcmfd.back();
}


void TransientSolver::setProlongation(bool prolong){
  _prolongation = prolong;
}








