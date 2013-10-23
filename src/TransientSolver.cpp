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
 
  /* set boolean paramters */
  _implicit = false;
  _improved = false;
  _steady_state = false;
  _adj_weight = false;

  /* compute the volume of the core */
  computeVolCore();

  /* compute new flux arrays */
  _mesh->createNewFlux(ADJ);
  _mesh->createNewFlux(PREVIOUS_CONV);
  _mesh->createNewFlux(FORWARD);
  _mesh->createNewFlux(REFERENCE);

}

/**
 * Transient Destructor clears all memory
 */
TransientSolver::~TransientSolver() {
}


void TransientSolver::solveInitialState(){

  log_printf(NORMAL, "Solving time dependent problem...");
  
  /* print problem specifications to screen */
  printProblemSpecs();
  
  /* initialize the time stepper and transient material properties */
  initializeTimeStepper();
  initializeTransientMaterials();

  /* compute the diffusion adjoint flux if adjoint weighting is turn on */
  if (_adj_weight){
    _k_eff_0 = _cmfd->computeKeff(ADJOINT);
    _mesh->copyFlux(CURRENT, ADJ);
    log_printf(NORMAL, "Adjoint k_eff: %f", _k_eff_0);
  }
  
  /* compute initial shape function */
  if (_solve_method == DIFFUSION)
    _k_eff_0 = _cmfd->computeKeff();
  else
    _k_eff_0 = _solver->convergeSource(1000);

  /* set the initial eigenvalue in the cmfd solver */
  log_printf(NORMAL, "k_eff_0: %f", _k_eff_0);
  _cmfd->setKeff0(_k_eff_0);
  
  /* store the initial flux and material properties */
  _mesh->copyFlux(CURRENT, PREVIOUS);
  _mesh->copyFlux(CURRENT, PREVIOUS_CONV);
  
  /* initialize the precursor concentration and PKE matrices */
  initializePrecursorConc();
  copyPrecConc(CURRENT, PREVIOUS);
  copyPrecConc(CURRENT, PREVIOUS_CONV);

  /* print out the initial N matrix */
  constructN();  
  for (int e = 0; e < _num_pke_groups; e++){
    log_printf(INFO, "e: %i, B: %f, C: %f", e, _B[e], _C[e]);
    for (int g = 0; g < _num_pke_groups; g++)
      log_printf(INFO, "g: %i, N: %f, Binter: %f", g, _N[e*_num_pke_groups+g], _B_inter[g]);
  }

  /* compute the power normalization factor */
  _power_factor = _power_init / computePower();
  
  /* save initial global parameters */
  _time.push_back(0.0);
  _temp.push_back(computeCoreTemp());
  _power.push_back(computePower() * _power_factor);
  log_printf(NORMAL, "TIME: %f, POWER: %.10f, TEMP: %f", 0.0, _power.back(), _temp.back());  
  _cmfd->setInitialSolve(false);
}


void TransientSolver::solveOuterStep(){
  
  /* set time integration parameters */
  double tolerance_out = 1.e-6;
  int length_conv = _power.size();
  double residual_out = 1.0;
  double keff;

  if (_ts->getTime(CURRENT) == _start_time)
    _timer->startTimer();
  
  /* compute the flux at the forward time step */
  if (_improved || _implicit){
    solvePKEsWithFeedback();
    updatePrecursorConc(PREVIOUS);
    keff = computeNewShape(FORWARD);
    copyFieldVariables(PREVIOUS_CONV, CURRENT);
    copyAmplitude(_B_prev, _B);
    _ts->setTime(CURRENT, _ts->getTime(PREVIOUS_CONV));
    _ts->setTime(PREVIOUS, _ts->getTime(PREVIOUS_CONV));
  }

  /* reset outer residual value and save current length of global vectors*/
  residual_out = 1.0;

  /* OUTER LOOP */
  while(residual_out > tolerance_out){
    
    trimVectors(length_conv);
    
    /* reset material properties and save old forward time step flux */
    if (_implicit == true && residual_out != 1.0){
      copyFieldVariables(PREVIOUS_CONV, CURRENT);
      copyAmplitude(_B_prev,_B);
      _ts->setTime(CURRENT, _ts->getTime(PREVIOUS_CONV));
      _ts->setTime(PREVIOUS, _ts->getTime(PREVIOUS_CONV));
    }

    /* solve PKE equations while perturbing materials and performing thermal feedback */
    solvePKEsWithFeedback();

    /* compute new forward time step shape */
    if (_implicit)
      updatePrecursorConc(FORWARD);
    else
      updatePrecursorConc(PREVIOUS);

    keff = computeNewShape();
    
    /* check the source residual */
    residual_out = computeResidual();
    if (_implicit)
      log_printf(NORMAL, "residual: %f", residual_out);
    
    copyFieldVariables(CURRENT, FORWARD);
  }

  _ts->setTime(PREVIOUS_CONV, _ts->getTime(CURRENT));
  _ts->incrementTime(FORWARD, _dt_outer);

  /* copy current solution vector to previous */
  copyAmplitude(_B,_B_prev);
  copyFieldVariables(CURRENT, PREVIOUS_CONV);
  copyFieldVariables(CURRENT, PREVIOUS);
  
  if (fabs(_ts->getTime(CURRENT) - _end_time) < 1e-4){
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


void TransientSolver::solvePKEsWithFeedback(){

  while(_ts->getTime(CURRENT) < _ts->getTime(FORWARD) - 1e-8){
    
    /* update materials and flux shape for current intermediate time step*/
    updateFlux();
    
    /* Compute PKE parameters */
    constructN();
    
    /* solver the PKE equations */
    solvePKEs();
      
    _ts->setTime(PREVIOUS, _ts->getTime(CURRENT));
    updateTemperatures();
       
    _time.push_back(_ts->getTime(CURRENT));
    _temp.push_back(computeCoreTemp());
    _power.push_back(computePower() * _power_factor);
    log_printf(NORMAL, "TIME: %f, POWER: %.10f, TEMP: %f", _time.back(), _power.back(), _temp.back());
  }
}

void TransientSolver::solvePKEs(){

  while(_ts->getTime(CURRENT) < _ts->getTime(PREVIOUS) + _dt_intermediate - 1e-9){
    
    copyAmplitude(_B, _B_inter);
    for (int e = 0; e < _num_pke_groups; e++){
      _B[e] += _C[e] * _dt_pke;
      for (int g = 0; g < _num_pke_groups; g++){
	_B[e] += _N[e*_num_pke_groups+g] * _B_inter[g] * _dt_pke;
      }
    }
    
    _ts->incrementTime(CURRENT, _dt_pke);
  }
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
  


void TransientSolver::copyAmplitude(double* b_from, double* b_to){

  for (int e = 0; e < _num_pke_groups; e++)
    b_to[e] = b_from[e];
}


/* compute new current flux by interpolating between previous_conv and forward flux */
void TransientSolver::updateFlux(){
  
  if (_improved){
    
    double ratio = _ts->getImprovedRatio();
    for (int i = 0; i < _mesh->getNumCells(); i++){
      for (int g = 0; g < _num_groups; g++){
	_mesh->getFluxes(CURRENT)[i*_num_groups + g] = (1 - ratio) * _mesh->getFluxes(PREVIOUS_CONV)[i*_num_groups + g]
	  + ratio * _mesh->getFluxes(FORWARD)[i*_num_groups + g];
	}
    }    
  }
  else if (_implicit)
    _mesh->copyFlux(FORWARD, CURRENT);

  normalizeFlux(CURRENT);
}


double TransientSolver::computeNewShape(materialState state){
  
  double keff = 0;
  
  syncMaterials(state);

  if (state != CURRENT)
    _mesh->copyFlux(CURRENT, REFERENCE);
  
  /* get keff */
  constructN();
  
  /* compute omega */
  double omega[_num_pke_groups];
  
  for (int e = 0; e < _num_pke_groups; e++)
    omega[e] = 0.0;
  
  for (int e = 0; e < _num_pke_groups; e++){
    omega[e] += _C[e];
    for (int g = 0; g < _num_pke_groups; g++)
      omega[e] += _N[e*_num_pke_groups+g] * _B[g];
  }
  
  for (int g = 0; g < _num_pke_groups; g++){
    _cmfd->setAmplitude(g, _B[g]);
    _cmfd->setFrequency(g, omega[g]);
  }
  
  if (_solve_method == DIFFUSION)
    keff = _cmfd->computeKeff();
  else
    keff = _solver->convergeSource(1000);
  
  /* copy flux to input state */
  if (state != CURRENT){
    _mesh->copyFlux(CURRENT, state);
    _mesh->copyFlux(REFERENCE, CURRENT);
  }

  log_printf(INFO, "time: %f, keff: %f", _ts->getTime(state), keff);
  
  return keff;
}


void TransientSolver::trimVectors(int len_conv){

  int len_curr = _power.size();
  for (int i = len_conv; i < len_curr; i++){
    _time.pop_back();
    _power.pop_back();
    _temp.pop_back();
  }
}

void TransientSolver::normalizeFlux(materialState state){
  
  double* flux = _mesh->getFluxes(state);
  
  /* normalize flux */
  if (_power_factor > 0.0){
    double power_old = _power.back();
    double power_new = computePower(state) * _power_factor;
    
    for (int i = 0; i < _mesh->getNumCells(); i++){
      for (int e = 0; e < _num_groups; e++){
	flux[i*_num_groups + e] = flux[i*_num_groups + e] * power_old / power_new;
      }
    }
  }  
}


double TransientSolver::computeResidual(){

  double res = 0.0;

  if (_implicit){
    double* new_flux = _mesh->getFluxes(CURRENT);
    double* fwd_flux = _mesh->getFluxes(FORWARD);
    
    for (int i = 0; i < _mesh->getNumCells(); i++){
      for (int e = 0; e < _num_groups; e++)
	res += pow(_materials[i]->getNuSigmaF()[e] * (fwd_flux[e] - new_flux[e]) * _mesh->getVolumes()[i], 2.0);
    }

    res = pow(res, 0.5);
    res = res / (_mesh->getNumCells() * _num_groups);
    
  }
  else
    res = 1e-10;

  return res;
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
    if (_multigroup_PKE){
      for (int e = 0; e < _num_groups; e++){
	power += _alpha / _nu * _materials[i]->getNuSigmaF()[e] * flux[i*_num_groups + e] * _B[e];
      }
    }
    else{
      for (int e = 0; e < _num_groups; e++){
	power += _alpha / _nu * _materials[i]->getNuSigmaF()[e] * flux[i*_num_groups + e] * _B[0];
      }
    }

    power = power * _power_factor;
    temp = _materials[i]->getTemperature(CURRENT) + _dt_intermediate * power;
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
  
  double fis_rate = 0.0;
  double global_conc[_num_delay_groups];
  
  for (int dg = 0; dg < _num_delay_groups; dg++){
    global_conc[dg] = 0.0;
  }
  
  double* flux = _mesh->getFluxes(CURRENT);
  
  for (int i = 0; i < _mesh->getNumCells(); i++){
    
    for (int dg = 0; dg < _num_delay_groups; dg++){
      fis_rate = 0.0;
      
      for (int e = 0; e < _num_groups; e++){
	fis_rate += _materials[i]->getNuSigmaF()[e] * flux[i*_num_groups + e] / _k_eff_0;
      }
    
      if (fis_rate > 0.0){
	static_cast<FunctionalMaterial*>(_materials[i])->setPrecConc(REFERENCE, _cmfd->getBeta()[dg] / _cmfd->getLambda()[dg] * fis_rate, dg);
	static_cast<FunctionalMaterial*>(_materials[i])->setPrecConc(PREVIOUS, _cmfd->getBeta()[dg] / _cmfd->getLambda()[dg] * fis_rate, dg);
	static_cast<FunctionalMaterial*>(_materials[i])->setPrecConc(PREVIOUS_CONV, _cmfd->getBeta()[dg] / _cmfd->getLambda()[dg] * fis_rate, dg);
	static_cast<FunctionalMaterial*>(_materials[i])->setPrecConc(CURRENT, _cmfd->getBeta()[dg] / _cmfd->getLambda()[dg] * fis_rate, dg);
	static_cast<FunctionalMaterial*>(_materials[i])->setPrecConc(FORWARD, _cmfd->getBeta()[dg] / _cmfd->getLambda()[dg] * fis_rate, dg);
	global_conc[dg] += _cmfd->getBeta()[dg] / _cmfd->getLambda()[dg] * fis_rate;
	log_printf(DEBUG, "cell: %i, dg: %i, initial prec: %f", i, dg, _cmfd->getBeta()[dg] / _cmfd->getLambda()[dg] * fis_rate);
      }
    }
  }
  
  for (int dg = 0; dg < _num_delay_groups; dg++){
    log_printf(INFO, "initial precursor conc: %f", global_conc[dg]);
  }
}


void TransientSolver::updatePrecursorConc(materialState state){

  log_printf(INFO, "Updating precursor concentrations...");
  
  /* initialize variables */
  double power, power_prev, old_conc, new_conc;
  double prefactor[3*_num_delay_groups];
  double global_conc[_num_delay_groups];
  
  //  double* flux = _mesh->getFluxes(CURRENT);
  double* flux = _mesh->getFluxes(state);
  
  /* compute prefactors */
  for (int i = 0; i < 3; i++){
    for (int dg = 0; dg < _num_delay_groups; dg++){
      global_conc[dg] = 0.0;
      prefactor[3*dg + 0] = exp(- _cmfd->getLambda()[dg] * _dt_outer);
      prefactor[3*dg + 1] = 1.0 - (1.0 - exp(- _cmfd->getLambda()[dg] * _dt_outer)) / (_cmfd->getLambda()[dg] * _dt_outer);
      prefactor[3*dg + 2] = exp(- _cmfd->getLambda()[dg] * _dt_outer) - (1.0 - exp(- _cmfd->getLambda()[dg] * _dt_outer)) / (_cmfd->getLambda()[dg] * _dt_outer);
    }
  }
  
  /* loop over fsrs */
  for (int i = 0; i < _mesh->getNumCells(); i++){
    
    for (int dg = 0; dg < _num_delay_groups; dg++){
      
      power = 0.0;
      power_prev = 0.0;
      
      /* compute power in fsr */
      if (_multigroup_PKE){
	for (int e = 0; e < _num_groups; e++){
	  power_prev += _materials[i]->getNuSigmaF()[e]*flux[i*_num_groups+e]*_B_prev[e] / _k_eff_0;
	  //power      += _materials[i]->getNuSigmaF()[e]*    flux[i*_num_groups+e]*_B[e]      / _k_eff_0;
	}
      }
      else{
	for (int e = 0; e < _num_groups; e++){
	  power_prev += _materials[i]->getNuSigmaF()[e]*flux[i*_num_groups+e]*_B_prev[0] / _k_eff_0;
	  //power      += _materials[i]->getNuSigmaF()[e]*    flux[i*_num_groups+e]*_B[0]      / _k_eff_0;
	}
      }
      
      /* compute new precursor concentration */
      //			old_conc = fsr_prev->getPrecursorConc()[dg];
      //			new_conc = prefactor[3*dg + 0]*old_conc + prefactor[3*dg + 1]*_beta[dg]/_lambda[dg]*power - prefactor[3*dg + 2]*_beta[dg]/_lambda[dg]*power_prev;
      //			fsr->setPrecursorConc(dg, new_conc);
      //			global_conc[dg] += new_conc;
      
      if (_materials[i]->isFissionable()){
	old_conc = static_cast<FunctionalMaterial*>(_materials[i])->getPrecConc(CURRENT, dg);
	new_conc = old_conc / (1 + _cmfd->getLambda()[dg]*_dt_outer) + _cmfd->getBeta()[dg]*_dt_outer*power_prev/(1.0+_cmfd->getLambda()[dg]*_dt_outer);
	log_printf(DEBUG, "cell: %i, old conc: %f, new conc: %f", i, old_conc, new_conc);
	static_cast<FunctionalMaterial*>(_materials[i])->setPrecConc(CURRENT, new_conc, dg);
	global_conc[dg] += new_conc;
      }
    }
  }
}


/* compute the power / cc */
double TransientSolver::computePower(materialState state){

  log_printf(INFO,"Computing power...");
  
  double power = 0.0;
  syncMaterials(state);
  double* flux = _mesh->getFluxes(state);
  double* volumes = _mesh->getVolumes();
  
  for (int i = 0; i < _mesh->getNumCells(); i++){
    
    if (_multigroup_PKE){
      for (int e = 0; e < _num_groups; e++){
	power += _kappa / _nu * _materials[i]->getNuSigmaF()[e]*flux[i*_num_groups + e] * volumes[i] * _B[e];
      }
    }
    else{
      for (int e = 0; e < _num_groups; e++){
	power += _kappa / _nu * _materials[i]->getNuSigmaF()[e]*flux[i*_num_groups + e] * volumes[i] * _B[0];
      }
    }
  }

  log_printf(DEBUG, "core power: %f, core power density: %f", power, power / _vol_core);
  
  power = power / _vol_core;
  
  return power;
}


void TransientSolver::constructN(){
  
  log_printf(INFO, "Constructing N...");
  
  syncMaterials(CURRENT);
  
  /* initialized variables */
  int indice1, indice2;
  double value;
  double beta_sum = 0.0;
  int ng = _num_groups;
  if (!_multigroup_PKE)
    ng = 1;
  
  double ratio = _ts->getImprovedRatio();
  
  /* get access to flux data */
  double* adjoint_flux;
  double* prev_flux = _mesh->getFluxes(PREVIOUS_CONV);
  double* mesh_flux = _mesh->getFluxes(CURRENT);
  if (_adj_weight)
    adjoint_flux = _mesh->getFluxes(ADJ);
  
  /* get access to mesh volumes */
  double* volumes = _mesh->getVolumes();
  
  /* get access to transient specific data */
  double* beta = _cmfd->getBeta();
  double* lambda = _cmfd->getLambda();
  double* velocities = _cmfd->getVelocity();
  
  int cw = _mesh->getCellsX();
  int ch = _mesh->getCellsY();
  
  double norm_factor[ng];
  int cell;
  
  /* zero normalization factor and leakage */
  for (int i = 0; i < ng; i++){
    norm_factor[i] = 0.0;
  }
  
  /* zero N matrix and C array */
  for (int e = 0; e < ng; e++){
    _C[e] = 0.0;
    for (int g = 0; g < ng; g++){
      _N[e*ng+g] = 0.0;
    }
  }
  
  /* compute the norm factor */
  for (int i = 0; i < _mesh->getNumCells(); i++){
    for (int e = 0; e < _num_groups; e++){
      if (_multigroup_PKE){
	if (_adj_weight)
	  norm_factor[e] += (volumes[i] * mesh_flux[i*_num_groups + e] * adjoint_flux[i*_num_groups + e]) / velocities[e];
	else
	  norm_factor[e] += (volumes[i] * mesh_flux[i*_num_groups + e]) / velocities[e];
      }
      else {
	if (_adj_weight)
	  norm_factor[0] += (volumes[i] * mesh_flux[i*_num_groups + e] * adjoint_flux[i*_num_groups + e]) / velocities[e];
	else
	  norm_factor[0] += (volumes[i] * mesh_flux[i*_num_groups + e]) / velocities[e];
      }
    }
  }
  
  /* take the inverse of the normalization factor */
  for (int e = 0; e < ng; e++){
    norm_factor[e] = 1.0 / norm_factor[e];
    
    log_printf(DEBUG, "norm factor group %i: %f", e, norm_factor[e]);
  }
  
  /* sum the beta */
  for (int dg = 0; dg < _num_delay_groups; dg++)
    beta_sum += beta[dg];

  log_printf(DEBUG, "beta sum: %f", beta_sum);
  
  /* loop over mesh cells */
  for (int y = 0; y < ch; y++){
    for (int x = 0; x < cw; x++){
      
      cell = y*cw+x;
      
      /* loop over energy groups */
      for (int e = 0; e < _num_groups; e++){
	
	/* shape derivative terms for IQS method */
	if (_transient_method == IQS){
	  if (ratio != 0){
	    if (_multigroup_PKE){
	      indice1 = e;
	      indice2 = e;
	      if (_adj_weight)
		value = 1.0/(velocities[e] * ratio * _dt_outer) * (prev_flux[cell*_num_groups + e] - mesh_flux[cell*_num_groups + e]) * volumes[cell] * norm_factor[e] * adjoint_flux[cell*_num_groups + e];
	      else
		value = 1.0/(velocities[e] * ratio * _dt_outer) * (prev_flux[cell*_num_groups + e] - mesh_flux[cell*_num_groups + e]) * volumes[cell] * norm_factor[e];
	    }
	    else{
	      indice1 = 0;
	      indice2 = 0;
	      if (_adj_weight)
		value = 1.0/(velocities[e] * ratio * _dt_outer) * (prev_flux[cell*_num_groups + e] - mesh_flux[cell*_num_groups + e]) * volumes[cell] * norm_factor[0] * adjoint_flux[cell*_num_groups + e];
	      else
		value = 1.0/(velocities[e] * ratio * _dt_outer) * (prev_flux[cell*_num_groups + e] - mesh_flux[cell*_num_groups + e]) * volumes[cell] * norm_factor[0];
	    }
	    
	    _N[indice1*_num_groups + indice2] += value;
	  }
	}
	/* absorption */
	if (_multigroup_PKE){
	  indice1 = e;
	  indice2 = e;
	  if (_adj_weight)
	    value = - _materials[cell]->getSigmaA()[e]*mesh_flux[cell*_num_groups + e]*volumes[cell] * norm_factor[e] * adjoint_flux[cell*_num_groups + e];
	  else
	    value = - _materials[cell]->getSigmaA()[e]*mesh_flux[cell*_num_groups + e]*volumes[cell] * norm_factor[e];
	}
	else{
	  indice1 = 0;
	  indice2 = 0;
	  if (_adj_weight)
	    value = - _materials[cell]->getSigmaA()[e]*mesh_flux[cell*_num_groups + e]*volumes[cell] * norm_factor[0] * adjoint_flux[cell*_num_groups + e];
	  else
	    value = - _materials[cell]->getSigmaA()[e]*mesh_flux[cell*_num_groups + e]*volumes[cell] * norm_factor[0];
	}
	
	_N[indice1*_num_groups + indice2] += value;
	
	/* loop over energy groups */
	for (int g = 0; g < _num_groups; g++){
	  
	  /* prompt fission for amplitude functions */
	  if (_multigroup_PKE){
	    indice1 = e;
	    indice2 = g;
	    if (_adj_weight)
	      value = (1-beta_sum)*_materials[cell]->getChi()[e]*_materials[cell]->getNuSigmaF()[g]*mesh_flux[cell*_num_groups + g]*volumes[cell]/_k_eff_0 * norm_factor[e] * adjoint_flux[cell*_num_groups + e];
	    else
	      value = (1-beta_sum)*_materials[cell]->getChi()[e]*_materials[cell]->getNuSigmaF()[g]*mesh_flux[cell*_num_groups + g]*volumes[cell]/_k_eff_0 * norm_factor[e];
	  }
	  else{
	    indice1 = 0;
	    indice2 = 0;
	    if (_adj_weight)
	      value = (1-beta_sum)*_materials[cell]->getChi()[e]*_materials[cell]->getNuSigmaF()[g]*mesh_flux[cell*_num_groups + g]*volumes[cell]/_k_eff_0 * norm_factor[0] * adjoint_flux[cell*_num_groups + e];
	    else
	      value = (1-beta_sum)*_materials[cell]->getChi()[e]*_materials[cell]->getNuSigmaF()[g]*mesh_flux[cell*_num_groups + g]*volumes[cell]/_k_eff_0 * norm_factor[0];
	  }
	  
	  _N[indice1*_num_groups + indice2] += value;

	  if (e != g){
	    /* out scattering */
	    if (_multigroup_PKE){
	      indice1 = e;
	      indice2 = e;
	      if (_adj_weight)
		value = - _materials[cell]->getSigmaS()[g*_num_groups + e]*mesh_flux[cell*_num_groups + e]*volumes[cell] * norm_factor[e] * adjoint_flux[cell*_num_groups + e];
	      else
		value = - _materials[cell]->getSigmaS()[g*_num_groups + e]*mesh_flux[cell*_num_groups + e]*volumes[cell] * norm_factor[e];
	    }
	    else{
	      indice1 = 0;
	      indice2 = 0;
	      if (_adj_weight)
		value = - _materials[cell]->getSigmaS()[g*_num_groups + e]*mesh_flux[cell*_num_groups + e]*volumes[cell] * norm_factor[0] * adjoint_flux[cell*_num_groups + e];
	      else
		value = - _materials[cell]->getSigmaS()[g*_num_groups + e]*mesh_flux[cell*_num_groups + e]*volumes[cell] * norm_factor[0];
	    }
	    
	    _N[indice1*_num_groups + indice2] += value;
	    
	    /* in scattering */
	    if (_multigroup_PKE){
	      indice1 = e;
	      indice2 = g;
	      if (_adj_weight)
		value = _materials[cell]->getSigmaS()[e*_num_groups + g]*mesh_flux[cell*_num_groups + g]*volumes[cell] * norm_factor[e] * adjoint_flux[cell*_num_groups + e];
	      else
		value = _materials[cell]->getSigmaS()[e*_num_groups + g]*mesh_flux[cell*_num_groups + g]*volumes[cell] * norm_factor[e];
	    }
	    else{
	      indice1 = 0;
	      indice2 = 0;
	      if (_adj_weight)
		value = _materials[cell]->getSigmaS()[e*_num_groups + g]*mesh_flux[cell*_num_groups + g]*volumes[cell] * norm_factor[0] * adjoint_flux[cell*_num_groups + e];
	      else
		value = _materials[cell]->getSigmaS()[e*_num_groups + g]*mesh_flux[cell*_num_groups + g]*volumes[cell] * norm_factor[0];
	    }
	    
	    _N[indice1*_num_groups + indice2] += value;
	  }
	}

	/* backwards difference fission for delayed precursors */
	for (int dg = 0; dg < _num_delay_groups; dg++){
	  for (int g = 0; g < _num_groups; g++){
	    if (_multigroup_PKE){
	      indice1 = e;
	      indice2 = g;
 	      if (_adj_weight)
		value = _dt_pke*lambda[dg]*beta[dg]*_materials[cell]->getChi()[e]*_materials[cell]->getNuSigmaF()[g]*mesh_flux[cell*_num_groups + g]*volumes[cell]/(_k_eff_0 * (1.0 + lambda[dg]*_dt_pke)) * norm_factor[e] * adjoint_flux[cell*_num_groups + e];
	      else
		value = _dt_pke*lambda[dg]*beta[dg]*_materials[cell]->getChi()[e]*_materials[cell]->getNuSigmaF()[g]*mesh_flux[cell*_num_groups + g]*volumes[cell]/(_k_eff_0 * (1.0 + lambda[dg]*_dt_pke)) * norm_factor[e];
	    }
	    else{
	      indice1 = 0;
	      indice2 = 0;
	      if (_adj_weight)
		value = _dt_pke*lambda[dg]*beta[dg]*_materials[cell]->getChi()[e]*_materials[cell]->getNuSigmaF()[g]*mesh_flux[cell*_num_groups + g]*volumes[cell]/(_k_eff_0 * (1.0 + lambda[dg]*_dt_pke)) * norm_factor[0] * adjoint_flux[cell*_num_groups + e];
	      else
		value = _dt_pke*lambda[dg]*beta[dg]*_materials[cell]->getChi()[e]*_materials[cell]->getNuSigmaF()[g]*mesh_flux[cell*_num_groups + g]*volumes[cell]/(_k_eff_0 * (1.0 + lambda[dg]*_dt_pke)) * norm_factor[0];
	      
	    }

	    _N[indice1*ng + indice2] += value;
	  }
	}
	
	/* previous time step precursor concentrations */
	if (_materials[cell]->getType() == FUNCTIONAL){
	  for (int dg = 0; dg < _num_delay_groups; dg++){
	    if (_multigroup_PKE){
	      indice1 = e;
	      if (_adj_weight)
		value = lambda[dg]*_materials[cell]->getChi()[e]*static_cast<FunctionalMaterial*>(_materials[cell])->getPrecConc(PREVIOUS, dg)*volumes[cell] / (1.0 + lambda[dg]*_dt_pke) * norm_factor[e] * adjoint_flux[cell*_num_groups + e];
	      else
		value = lambda[dg]*_materials[cell]->getChi()[e]*static_cast<FunctionalMaterial*>(_materials[cell])->getPrecConc(PREVIOUS, dg)*volumes[cell] / (1.0 + lambda[dg]*_dt_pke) * norm_factor[e];
	    }
	    else{
	      indice1 = 0;
	      if (_adj_weight)
		value = lambda[dg]*_materials[cell]->getChi()[e]*static_cast<FunctionalMaterial*>(_materials[cell])->getPrecConc(PREVIOUS, dg)*volumes[cell] / (1.0 + lambda[dg]*_dt_pke) * norm_factor[0] * adjoint_flux[cell*_num_groups + e];
	      else
		value = lambda[dg]*_materials[cell]->getChi()[e]*static_cast<FunctionalMaterial*>(_materials[cell])->getPrecConc(PREVIOUS, dg)*volumes[cell] / (1.0 + lambda[dg]*_dt_pke) * norm_factor[0];
	    }
	  _C[indice1] += value;
	  }
	}
      }
    }
  }

  
  for (int e = 0; e < _num_groups; e++){
    /* set leakage terms */
    value = _mesh->getLeakage(CURRENT, e, _adj_weight);

    if (_multigroup_PKE){
      _N[e*_num_groups + e] += value * norm_factor[e];
      log_printf(INFO, "g: %i, net leakage: %f", e, value * norm_factor[e]);
    }
    else
      _N[e*_num_groups + e] += value * norm_factor[0];
  }
}


void TransientSolver::setStartTime(double time){
  _start_time = time;
}


void TransientSolver::setEndTime(double time){
  _end_time = time;
}


void TransientSolver::initializeTimeStepper(){
  _ts = new TimeStepper(_start_time, _end_time, _dt_outer, _dt_intermediate);
  _cmfd->initialize(_ts, _adj_weight);
}


void TransientSolver::initializeTransientMaterials(){

  for (int i = 0; i < _mesh->getNumCells(); i++){
    if (_materials[i]->getType() == FUNCTIONAL){
      static_cast<FunctionalMaterial*>(_materials[i])->setTimeStepper(_ts);
      static_cast<FunctionalMaterial*>(_materials[i])->initializeTransientProps(_num_delay_groups);
    }
  }
}


void TransientSolver::setAdjointWeighting(bool weighting){
  _adj_weight = weighting;
}


void TransientSolver::setDtOuter(double dt){
  _dt_outer = dt;
}


void TransientSolver::setDtIntermediate(double dt){
  _dt_intermediate = dt;
}


void TransientSolver::setDtPKE(double dt){
  _dt_pke = dt;
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
}

void TransientSolver::setTransientMethod(const char* trans_type){

  if (strcmp("ADIABATIC", trans_type) == 0)
    _transient_method = ADIABATIC;
  else if (strcmp("OMEGA_MODE", trans_type) == 0)
    _transient_method = OMEGA_MODE;
  else if (strcmp("IQS", trans_type) == 0)
    _transient_method = IQS;
  else
    log_printf(ERROR, "Could not recognize transient method; "
	       " the options are ADIABATIC, OMEGA_MODE, and IQS");

  _cmfd->setTransientType(_transient_method);
}

void TransientSolver::setOuterImplicit(bool implicit){
  _implicit = implicit;
}

void TransientSolver::setImproved(bool improved){
  _improved = improved;
}

void TransientSolver::setMultigroupPKE(bool multigroup){

  _multigroup_PKE = multigroup;

  if (multigroup)
    _num_pke_groups = _num_groups;
  else
    _num_pke_groups = 1;

  /* allocate memory for PKE matrices and arrays */
  _N = new double[_num_pke_groups*_num_pke_groups];
  _B = new double[_num_pke_groups];
  _B_prev = new double[_num_pke_groups];
  _B_inter = new double[_num_pke_groups];
  _C = new double[_num_pke_groups];
  
  /* initialize the amplitude for each group to 1.0 */
  for (int e = 0; e < _num_pke_groups; e++){
    _B[e] = 1.0;
    _B_prev[e] = 1.0;
    _B_inter[e] = 1.0;
  }
}

void TransientSolver::setSteadyState(bool steady_state){
  _steady_state = steady_state;
  _cmfd->setMultigroupPKE(steady_state);
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


/* Print the problem specifications */
void TransientSolver::printProblemSpecs(){
  
  if (_transient_method == IQS)
    log_printf(NORMAL, "Solving transient problem with IQS method");
  else if (_transient_method == OMEGA_MODE)
    log_printf(NORMAL, "Solving transient problem with OMEGA MODE method");
  else if (_transient_method == ADIABATIC)
    log_printf(NORMAL, "Solving transient problem with ADIABATIC method");
  
  if (_implicit)
    log_printf(NORMAL,"IMPLICIT OUT --- TRUE");
  else
    log_printf(NORMAL,"IMPLICIT OUT --- FALSE");
  
  if (_improved)
    log_printf(NORMAL,"IMPROVED     --- TRUE");
  else
    log_printf(NORMAL,"IMPROVED     --- FALSE");
  
  if (_steady_state)
    log_printf(NORMAL,"STEADY STATE --- TRUE");
  else
    log_printf(NORMAL,"STEADY STATE --- FALSE");
  
  if (_multigroup_PKE)
    log_printf(NORMAL,"MG PKE       --- TRUE");
  else
    log_printf(NORMAL,"MG PKE       --- FALSE");
  
  if (_adj_weight)
    log_printf(NORMAL,"ADJ WEIGHT   --- TRUE");
  else
    log_printf(NORMAL,"ADJ WEIGHT   --- FALSE");
}
