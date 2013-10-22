/*
 * Cmfd.cpp
 *
 *  Created on: September 8, 2013
 *      Author: Sam Shaner
 *				MIT, Course 22
 *              shaner@mit.edu
 */

#include "CmfdTransient.h"

/**
 * CmfdTransient constructor
 * @param geom pointer to the geometry
 * @param criteria convergence criterial on keff 
 */
CmfdTransient::CmfdTransient(Geometry* geometry, double criteria) :
  Cmfd(geometry, criteria){
  
  /* Boolean and Enum flags to toggle features */
  _transient_method = NONE;
  
  /* Create Cmfd matrix objects */
  int petsc_err;
  petsc_err = createBPrime();
  _num_groups = geometry->getNumEnergyGroups();
  
  /* transient parameters */
  _lambda = NULL;
  _beta = NULL;
  _velocity = NULL;
  _frequency = NULL;
  _amplitude = NULL;
  _initial_solve = true;
  _k_eff_0 = 1.0;
  
}

/**
 * cmfd Destructor clears all memory
 */
CmfdTransient::~CmfdTransient() {
}


int CmfdTransient::createBPrime(){

  int petsc_err = 0;
  
  /* Create flux, source, and residual vectors */
  petsc_err = VecCreateSeq(PETSC_COMM_WORLD, _cells_x*_cells_y*_num_groups, &_b_prime);
  petsc_err = VecCreateSeq(PETSC_COMM_WORLD, _cells_x*_cells_y*_num_groups, &_b);
  CHKERRQ(petsc_err);
  
  return petsc_err;
}


/*
 * CMFD solver that solves the diffusion problem
 * @param solve methed - either diffusion or cmfd (acceleration)
 * @param iteration number of in MOC solver - used for plotting
 * @return k-effective
 */
double CmfdTransient::computeKeff(fluxType flux_method){

  log_printf(INFO, "Running cmfd diffusion diffusion solver...");
  
  _flux_method = flux_method;

  if (_solve_method == MOC)
    computeXS();

  computeDs();

  /* initialize variables */
  int petsc_err = 0;
  int iter = 0;
  PetscScalar sumold, sumnew, scale_val, eps;
  PetscReal rtol = 1e-10;
  PetscReal atol = 1e-10;
  Vec sold, snew, res;
  KSP ksp;
  
  /* zero matrices */
  petsc_err = MatZeroEntries(_A);
  petsc_err = VecZeroEntries(_b_prime);
  
  if (_assemble_M)
    petsc_err = MatZeroEntries(_M);

  if (_transient_method == IQS)
    petsc_err = VecZeroEntries(_b);
  
  /* construct matrices */
  petsc_err = constructMatrices();
  CHKERRQ(petsc_err);
    
  /* construct matrix and vector objects */
  petsc_err = VecAssemblyBegin(_phi_new);
  petsc_err = VecAssemblyEnd(_phi_new);
  petsc_err = VecAssemblyBegin(_phi_old);
  petsc_err = VecAssemblyEnd(_phi_old);
  petsc_err = VecAssemblyBegin(_b);
  petsc_err = VecAssemblyEnd(_b);
  petsc_err = VecAssemblyBegin(_sold);
  petsc_err = VecAssemblyEnd(_sold);
  petsc_err = VecAssemblyBegin(_snew);
  petsc_err = VecAssemblyEnd(_snew);
  petsc_err = VecAssemblyBegin(_res);
  petsc_err = VecAssemblyEnd(_res);
  petsc_err = VecAssemblyBegin(_b_prime);
  petsc_err = VecAssemblyEnd(_b_prime);
  petsc_err = MatAssemblyBegin(_A, MAT_FINAL_ASSEMBLY);
  petsc_err = MatAssemblyEnd(_A, MAT_FINAL_ASSEMBLY);

  /* assemble M in needed */
  if (_assemble_M){
    petsc_err = MatAssemblyBegin(_M, MAT_FINAL_ASSEMBLY);
    petsc_err = MatAssemblyEnd(_M, MAT_FINAL_ASSEMBLY);
  }
  CHKERRQ(petsc_err);

  /* create petsc ksp objects */
  petsc_err = KSPCreate(PETSC_COMM_WORLD, &ksp);
  petsc_err = KSPSetTolerances(ksp, rtol, atol, PETSC_DEFAULT, PETSC_DEFAULT);
  petsc_err = KSPSetType(ksp, KSPGMRES);
  petsc_err = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
  petsc_err = KSPSetOperators(ksp, _A, _A, SAME_NONZERO_PATTERN);
  petsc_err = KSPSetUp(ksp);
  petsc_err = KSPSetFromOptions(ksp);
  CHKERRQ(petsc_err);  
  
  /* If not solving a transport transient with IQS, use SLEPc eigenvalue solver */
  if (_transient_method != IQS){
    
    /* compute the normalize the initial source */
    petsc_err = MatMult(_M, _phi_old, _sold);
    petsc_err = VecSum(_sold, &sumold);
    scale_val = (_cells_x * _cells_y * _num_groups) / sumold;
    petsc_err = VecScale(_sold, scale_val);
    sumold = _cells_x * _cells_y * _num_groups;
    CHKERRQ(petsc_err);
    
    /* power iteration diffusion solver */
    for (iter = 0; iter < 1000; iter++){
      
      /* Solve phi = A^-1 * old_source */
      petsc_err = KSPSolve(ksp, _sold, _phi_new);
      
      /* compute the new source */
      petsc_err = MatMult(_M, _phi_new, _snew);
      petsc_err = VecSum(_snew, &sumnew);
      CHKERRQ(petsc_err);
      
      /* compute and set keff */
      _k_eff = sumnew / sumold;
      
      /* scale the old source by keff */
      petsc_err = VecScale(_sold, _k_eff);
      
      /* compute the L2 norm of source error */
      scale_val = 1e-15;
      petsc_err = VecShift(_snew, scale_val);
      petsc_err = VecShift(_sold, scale_val);
      petsc_err = VecPointwiseDivide(_res, _sold, _snew);
      scale_val = -1;
      petsc_err = VecShift(_res, scale_val);
      scale_val = -1e-15;
      petsc_err = VecShift(_snew, scale_val);
      petsc_err = VecShift(_sold, scale_val);
      CHKERRQ(petsc_err);
      petsc_err = VecNorm(_res, NORM_2, &eps);
      eps = eps / (_cells_x * _cells_y * _num_groups);
      
      /* normalize the new source */
      scale_val = (_cells_x * _cells_y * _num_groups) / sumnew;
      petsc_err = VecScale(_snew, scale_val);
      CHKERRQ(petsc_err);
      
      /* set old source to new source */
      petsc_err = VecCopy(_snew, _sold);
      CHKERRQ(petsc_err);
      
      log_printf(INFO, "CMFD iter: %i, keff: %f, error: %f", iter + 1, _k_eff, eps);
      
      /* check for convergence */
      if (eps < _conv_criteria)
	break;
    }
  }
  else{

    petsc_err = VecCopy(_phi_old, _phi_new);
    
    /* get initial source and find initial k_eff */
    petsc_err = MatMult(_M, _phi_new, _snew);
    petsc_err = VecSum(_snew, &sumnew);
    
    /* normalize b prime */
    scale_val = (_cells_x*_cells_y*_num_groups) / sumnew;
    petsc_err = VecScale(_b_prime, scale_val);
    CHKERRQ(petsc_err);
    
    /* compute A * phi + b_prime */
    petsc_err = VecScale(_b_prime, -1);
    petsc_err = MatMultAdd(_A, _phi_old, _b_prime, _b);
    petsc_err = VecScale(_b_prime, -1);
    petsc_err = VecSum(_b, &sumold);

    /* compute keff */
    _k_eff = float(sumnew)/float(sumold);
    log_printf(INFO, "CMFD iter: %i, keff: %.10f, snew: %f, sold: %f, scale_val: %f", iter, _k_eff, sumnew, sumold, scale_val);
    CHKERRQ(petsc_err);
    
    /* recompute and normalize initial source */
    petsc_err = MatMult(_M, _phi_old, _sold);
    petsc_err = VecSum(_sold, &sumold);
    scale_val = (_cells_x*_cells_y * _num_groups) / sumold;
    petsc_err = VecScale(_sold, scale_val);
    sumold = _cells_x*_cells_y * _num_groups;
    CHKERRQ(petsc_err);
    
    /* perform power iterations to converge the flux */
    for (iter = 0; iter < 1000; iter++){
      
      /* solve A * b = phi_new */
      petsc_err = VecWAXPY(_b, _k_eff, _b_prime, _sold);
      petsc_err = KSPSolve(ksp, _b, _phi_new);
      
      /* computed the new source */
      petsc_err = MatMult(_M, _phi_new, _snew);
      petsc_err = VecSum(_snew, &sumnew);
      
      /* compute new keff */
      _k_eff = sumnew / sumold;

      /* compute the L2 norm of source error */
      scale_val = 1e-15;
      petsc_err = VecShift(_snew, scale_val);
      petsc_err = VecShift(_sold, scale_val);
      petsc_err = VecPointwiseDivide(_res, _sold, _snew);
      scale_val = -1;
      petsc_err = VecShift(_res, scale_val);
      scale_val = -1e-15;
      petsc_err = VecShift(_snew, scale_val);
      petsc_err = VecShift(_sold, scale_val);
      CHKERRQ(petsc_err);
      petsc_err = VecNorm(_res, NORM_2, &eps);
      eps = eps / (_cells_x * _cells_y * _num_groups);
      
      /* normalize the new source */
      scale_val = (_cells_x*_cells_y * _num_groups) / sumnew;
      petsc_err = VecScale(_snew, scale_val);
      CHKERRQ(petsc_err);
      
      /* set old source to new source */
      petsc_err = VecCopy(_snew, _sold);
      CHKERRQ(petsc_err);
      
      /* check for convergence */
      if (eps < _conv_criteria)
	break;   
    }
  }
  
  /* destroy KSP object */
  petsc_err = KSPDestroy(&ksp);
  
  /* rescale flux and pass to meshCells */
  petsc_err = rescaleFlux();
  
  /* give the petsc flux array to the mesh cell flux array */
  petsc_err = setMeshCellFlux();
  
  updateMOCFlux();
  
  return _k_eff;
}


/* rescale flux */
int CmfdTransient::rescaleFlux(){

  /* initialize variables */
  int petsc_err = 0;
  PetscScalar sumnew, sumold, scale_val;

  /* rescale the new and old flux to have an avg source of 1.0 */
  petsc_err = VecSum(_phi_new, &sumnew);
  scale_val = _cells_x*_cells_y*_num_groups / sumnew;
  petsc_err = VecScale(_phi_new, scale_val);

  CHKERRQ(petsc_err);
  
  return 0;
}


/* Fill in the values in the A matrix, M matrix, and phi_old vector
 * @param A matrix
 * @param M matrix
 * @param old flux vector
 * @param solve methed - either DIFFUSION or CMFD
 * @return petsc error indicator
 */
int CmfdTransient::constructMatrices(){

  log_printf(INFO,"Constructing cmfd transient matrices...");
  
  /* initialized variables */
  int petsc_err = 0;
  PetscInt row, col;
  PetscScalar value, phi, b_prime = 0;
  int cell;
  
  Material **materials = _mesh->getMaterials();
  double* old_flux = _mesh->getFluxes(PREVIOUS);
  double* heights = _mesh->getLengthsY();
  double* widths = _mesh->getLengthsX();
  
  /* loop over mesh cells in y direction */
  for (int y = 0; y < _cells_y; y++){
    
    /* loop over mesh cells in x direction */
    for (int x = 0; x < _cells_x; x++){
      
      cell = y*_cells_x + x;
      
      /* loop over energy groups */
      for (int e = 0; e < _num_groups; e++){

	row = cell*_num_groups+e;

	/* flux */
	value = old_flux[cell*_num_groups + e];
	petsc_err = VecSetValue(_phi_old, row, value, INSERT_VALUES);
	
	/* absorption term */
	col = cell*_num_groups+e;
	value = materials[cell]->getSigmaA()[e] * _mesh->getVolumes()[cell];
	petsc_err = MatSetValues(_A, 1, &row, 1, &col, &value, ADD_VALUES);

	/* out (1st) and in (2nd) scattering */
	if (_flux_method == PRIMAL){
	  for (int g = 0; g < _num_groups; g++){
	    if (e != g){
	      col = cell*_num_groups+e;
	      value = materials[cell]->getSigmaS()[g*_num_groups + e] * _mesh->getVolumes()[cell]; 
	      petsc_err = MatSetValues(_A, 1, &row, 1, &col, &value, ADD_VALUES);

	      if (_multigroup_PKE)
		value = - materials[cell]->getSigmaS()[e*_num_groups + g] * _mesh->getVolumes()[cell] * _amplitude[g] / _amplitude[e];
	      else
		value = - materials[cell]->getSigmaS()[e*_num_groups + g] * _mesh->getVolumes()[cell] * _amplitude[0] / _amplitude[0];  

	      col = cell*_num_groups+g;
	      petsc_err = MatSetValues(_A, 1, &row, 1, &col, &value, ADD_VALUES);
	    }
	  }
	}
	else{
	  for (int g = 0; g < _num_groups; g++){
	    if (e != g){

	      col = cell*_num_groups+e;
	      value = materials[cell]->getSigmaS()[g*_num_groups + e] * _mesh->getVolumes()[cell];
	      petsc_err = MatSetValues(_A, 1, &col, 1, &row, &value, ADD_VALUES);

	      col = cell*_num_groups+g;
	      value = - materials[cell]->getSigmaS()[e*_num_groups + g] * _mesh->getVolumes()[cell];
	      petsc_err = MatSetValues(_A, 1, &col, 1, &row, &value, ADD_VALUES);
	    }
	  }
	}

	/* add frequency term to diagonal of A */
	if (_transient_method == OMEGA_MODE){
	  if (_multigroup_PKE)
	    value = _frequency[e] / (_velocity[e] * _amplitude[e]) * _mesh->getVolumes()[cell];
	  else
	    value = _frequency[0] / (_velocity[e] * _amplitude[0]) * _mesh->getVolumes()[cell];

	  col = cell*_num_groups+e;
	  petsc_err = MatSetValues(_A, 1, &row, 1, &col, &value, ADD_VALUES);
	}
	else if (_transient_method == IQS && !_initial_solve){
	  if (_multigroup_PKE)
	    value = 1.0 / _velocity[e] * (_frequency[e] / _amplitude[e] + 1.0 / _time_stepper->getDtOuter()) * _mesh->getVolumes()[cell];
	  else
	    value = 1.0 / _velocity[e] * (_frequency[0] / _amplitude[0] + 1.0 / _time_stepper->getDtOuter()) * _mesh->getVolumes()[cell];
	  
	  col = cell*_num_groups+e;
	  petsc_err = MatSetValues(_A, 1, &row, 1, &col, &value, ADD_VALUES);
	}
	
	/* add previous shape to new vector */
	if (_transient_method == IQS && !_initial_solve){
	  value = 1.0 / (_velocity[e] * _time_stepper->getDtOuter()) * _mesh->getVolumes()[cell] * old_flux[cell*_num_groups + e];
	  petsc_err = VecSetValue(_b_prime, row, value, INSERT_VALUES);
	}
		
	/* RIGHT SURFACE */
       
	/* set transport term on diagonal */
	value = (materials[cell]->getDifHat()[2*_num_groups + e] 
			    - materials[cell]->getDifTilde()[2*_num_groups + e]) 
	                    * heights[cell / _cells_x];

	col = cell*_num_groups+e;
	petsc_err = MatSetValues(_A, 1, &row, 1, &col, &value, ADD_VALUES);
	
	/* set transport term on off diagonal */
	if (x != _cells_x - 1){
	  value = - (materials[cell]->getDifHat()[2*_num_groups + e] 
					+ materials[cell]->getDifTilde()[2*_num_groups + e]) 
	                                * heights[cell / _cells_x];

	  col = (cell+1)*_num_groups+e;
	  petsc_err = MatSetValues(_A, 1, &row, 1, &col, &value, ADD_VALUES);
	}

	/* LEFT SURFACE */
	
	/* set transport term on diagonal */
	value = (materials[cell]->getDifHat()[0*_num_groups + e] 
			    + materials[cell]->getDifTilde()[0*_num_groups + e]) 
	                    * heights[cell / _cells_x];
	
	col = cell*_num_groups+e;
	petsc_err = MatSetValues(_A, 1, &row, 1, &col, &value, ADD_VALUES);
	
	/* set transport term on off diagonal */
	if (x != 0){
	    value = - (materials[cell]->getDifHat()[0*_num_groups + e] 
			    - materials[cell]->getDifTilde()[0*_num_groups + e]) 
	                    * heights[cell / _cells_x];

	    col = (cell-1)*_num_groups+e;
	  petsc_err = MatSetValues(_A, 1, &row, 1, &col, &value, ADD_VALUES);
	}
	
	/* BOTTOM SURFACE */
	
	/* set transport term on diagonal */
	value = (materials[cell]->getDifHat()[1*_num_groups + e] 
			    - materials[cell]->getDifTilde()[1*_num_groups + e]) 
	                    * widths[cell % _cells_x];

	col = cell*_num_groups+e;
	petsc_err = MatSetValues(_A, 1, &row, 1, &col, &value, ADD_VALUES);
       
	/* set transport term on off diagonal */
	if (y != _cells_y - 1){
	  value = - (materials[cell]->getDifHat()[1*_num_groups + e] 
			  + materials[cell]->getDifTilde()[1*_num_groups + e]) 
	                  * widths[cell % _cells_x];

	  col = (cell+_cells_x)*_num_groups+e;
	  petsc_err = MatSetValues(_A, 1, &row, 1, &col, &value, ADD_VALUES);
	}
	
	/* TOP SURFACE */
	
	/* set transport term on diagonal */
        value = (materials[cell]->getDifHat()[3*_num_groups + e] 
			    + materials[cell]->getDifTilde()[3*_num_groups + e]) 
	                    * widths[cell % _cells_x];

	col = cell*_num_groups+e;
	petsc_err = MatSetValues(_A, 1, &row, 1, &col, &value, ADD_VALUES);
	
	/* set transport term on off diagonal */
	if (y != 0){
	  value = - (materials[cell]->getDifHat()[3*_num_groups + e] 
					- materials[cell]->getDifTilde()[3*_num_groups + e]) 
	                                * widths[cell % _cells_x];

	  col = (cell-_cells_x)*_num_groups+e;
	  petsc_err = MatSetValues(_A, 1, &row, 1, &col, &value, ADD_VALUES);
	}
			
	/* add fission terms to M */
	if (_assemble_M){
	  for (int g = 0; g < _num_groups; g++){

	    if (_transient_method == ADIABATIC){
	      if (_multigroup_PKE)
		value = materials[cell]->getChi()[e] * materials[cell]->getNuSigmaF()[g] * _mesh->getVolumes()[cell] * _amplitude[g] / _amplitude[e] / _k_eff_0;
	      else
		value = materials[cell]->getChi()[e] * materials[cell]->getNuSigmaF()[g] * _mesh->getVolumes()[cell] * _amplitude[0] / _amplitude[0] / _k_eff_0;
	    }
	    else{
	      value = 0.0;
	      for (int dg = 0; dg < _num_delay_groups; dg++){
		if (_multigroup_PKE)
		  value += (_beta[dg] * _lambda[dg] * materials[cell]->getChi()[e]* materials[cell]->getNuSigmaF()[g]) / ((_lambda[dg] + _frequency[dg]) * _k_eff_0) * _amplitude[g] / _amplitude[e];
		else
		  value += (_beta[dg] * _lambda[dg] * materials[cell]->getChi()[e]* materials[cell]->getNuSigmaF()[g]) / ((_lambda[dg] + _frequency[dg]) * _k_eff_0) * _amplitude[0] / _amplitude[0];
	      }

	      if (_multigroup_PKE)
		value += materials[cell]->getChi()[e] * materials[cell]->getNuSigmaF()[g] * _mesh->getVolumes()[cell] * (1.0 - _beta_sum)  * _amplitude[g] / _amplitude[e] / _k_eff_0;
	      else
		value += materials[cell]->getChi()[e] * materials[cell]->getNuSigmaF()[g] * _mesh->getVolumes()[cell] * (1.0 - _beta_sum)  * _amplitude[0] / _amplitude[0] / _k_eff_0;	      
	    }

	    col = cell*_num_groups+g;
	    if (_flux_method == PRIMAL)
	      petsc_err = MatSetValues(_M, 1, &row, 1, &col, &value, ADD_VALUES);
	    else
	      petsc_err = MatSetValues(_M, 1, &col, 1, &row, &value, ADD_VALUES);	    
	  }
	}
      }
    }
  }
  
  log_printf(INFO,"Done constructing matrices...");
  
  return petsc_err;
}


void CmfdTransient::setFrequency(int e, double value){
  _frequency[e] = value;
}

void CmfdTransient::setAmplitude(int e, double value){
  _amplitude[e] = value;
}

void CmfdTransient::setTransientType(transientType trans_type){
  _transient_method = trans_type;
}

transientType CmfdTransient::getTransientType(){
  return _transient_method;
}

void CmfdTransient::setKeff0(double keff_0){
  _k_eff_0 = keff_0;
}

void CmfdTransient::setMultigroupPKE(bool mg){

  _multigroup_PKE = mg;
  
  if (_frequency != NULL)
    delete [] _frequency;
  
  if (_amplitude != NULL)
    delete [] _amplitude;

  /* allocate memory and initialize frequency and amplitude arrays */
  if (mg){
    _frequency = new double[_num_groups];
    _amplitude = new double[_num_groups];
    for (int i = 0; i < _num_groups; i++){
      _frequency[i] = 0.0;
      _amplitude[i] = 1.0;
    }
  }
  else{
    _frequency = new double[1];
    _amplitude = new double[1];
    _frequency[0] = 0.0;
    _amplitude[0] = 1.0;
  }


}


void CmfdTransient::initialize(TimeStepper* time_stepper, bool adjoint_weighted){

  _adj_weight = adjoint_weighted;
  Material **materials = _mesh->getMaterials();
  _time_stepper = time_stepper;

  if (_solve_method == MOC){
    for (int i = 0; i < _num_fsrs; i++){
      if (_FSR_materials[i]->getType() == FUNCTIONAL){
      static_cast<FunctionalMaterial*>(_FSR_materials[i])->setTimeStepper(time_stepper);
      }
    }
  }
}


double* CmfdTransient::getBeta(){
  return _beta;
}

double* CmfdTransient::getLambda(){
  return _lambda;
}

double* CmfdTransient::getVelocity(){
  return _velocity;
}

void CmfdTransient::setBeta(double* beta, int num_delay_groups){

  _beta = new double[num_delay_groups];
  
  for (int g = 0; g < num_delay_groups; g++){
    _beta[g] = beta[g];
  }
}


void CmfdTransient::setLambda(double* decay_const, int num_delay_groups){

  _lambda = new double[num_delay_groups];
  
  for (int g = 0; g < num_delay_groups; g++){
    _lambda[g] = decay_const[g];
  }
  
}


void CmfdTransient::setVelocity(double* velocity, int num_groups){

  _velocity = new double[num_groups];
  
  for (int g = 0; g < num_groups; g++)
    _velocity[g] = velocity[g];
  
}


int CmfdTransient::dumpVec(Vec petsc_vec){

  log_printf(NORMAL, "Dumping vector...");

  int petsc_err = 0;

  PetscScalar *scalar_vec;
  petsc_err = VecGetArray(petsc_vec, &scalar_vec);
  CHKERRQ(petsc_err);
  
  for (int i = 0; i < _cells_x*_cells_y; i++){
    for (int e = 0; e < _num_groups; e++){
      log_printf(NORMAL, "cell: %i, group: %i, value: %f", i, e, double(scalar_vec[i*_num_groups + e]));
    }
  }
  
  petsc_err = VecRestoreArray(petsc_vec, &scalar_vec);
  CHKERRQ(petsc_err);

  return petsc_err;
}


void CmfdTransient::setInitialSolve(bool initial){
  _initial_solve = initial;
}
