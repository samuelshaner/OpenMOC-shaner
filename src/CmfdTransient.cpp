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
  _omega = 1.0;

  /* set number of groups */
  _num_groups = geometry->getNumEnergyGroups();

  /* create vector objects */
  _b_prime = new double[_cells_x*_cells_y*_num_groups];
  _b = new double[_cells_x*_cells_y*_num_groups];
  
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

  delete [] _b;
  delete [] _b_prime;

}



/*
 * CMFD solver that solves the diffusion problem
 * @param solve methed - either diffusion or cmfd (acceleration)
 * @param iteration number of in MOC solver - used for plotting
 * @return k-effective
 */
double CmfdTransient::computeKeff(fluxType flux_method){

  log_printf(INFO, "Running cmfd diffusion diffusion solver...");
  
  /* initialize variables */
  int iter = 0;
  double norm = 0.0;
  double sumold, sumnew, scale_val, eps;
  double val = 0.0;
  int row = 0;
  int iter_in = 0;
  double conv = 1e1;

  _flux_method = flux_method;

  if (_solve_method == MOC)
    computeXS();

  computeDs();

  /* construct matrices */
  constructMatrices();

  vecCopy(_phi_old, _phi_new);
  vecCopy(_phi_old, _phi_temp);
    
  /* If not solving a transport transient with IQS, use SLEPc eigenvalue solver */
  if (_transient_method != IQS){
    
    /* compute the normalize the initial source */
    matMult(_M, _phi_new, _sold);
    sumold = vecSum(_sold);
    scale_val = (_cells_x * _cells_y * _num_groups) / sumold;
    vecScale(_sold, scale_val);
    sumold = _cells_x * _cells_y * _num_groups;
    
    /* power iteration diffusion solver */
    for (iter = 0; iter < 20000; iter++){
      
      /* solver phi = A^-1 * old_source */
      linearSolve(_A, _phi_new, _sold, conv, _omega);

      /* compute the new source */
      matMult(_M, _phi_new, _snew);
      sumnew = vecSum(_snew);
      
      /* compute and set keff */
      _k_eff = sumnew / sumold;
      
      /* scale the old source by keff */
      vecScale(_sold, _k_eff);
      
      /* compute the L2 norm of source error */
      norm = 0.0;
      for (int i = 0; i < _cells_x*_cells_y*_num_groups; i++)
	norm += pow(_snew[i] - _sold[i], 2);
      
      norm = pow(norm, 0.5);
      norm = norm / (_cells_x*_cells_y*_num_groups);
      scale_val = (_cells_x*_cells_y*_num_groups) / sumnew;
      vecScale(_snew, scale_val);
      vecCopy(_snew, _sold);
      
      /* set old source to new source */
      vecCopy(_snew, _sold);
      
      log_printf(INFO, "CMFD iter: %i, keff: %f, error: %f", iter + 1, _k_eff, norm);
      
      /* check for convergence */
      if (norm < _conv_criteria)
	break;
    }
  }
  else{

    vecCopy(_phi_old, _phi_new);
    
    /* get initial source and find initial k_eff */
    matMult(_M, _phi_new, _snew);
    sumnew = vecSum(_snew);
    
    /* normalize b prime */
    scale_val = (_cells_x*_cells_y*_num_groups) / sumnew;
    vecScale(_b_prime, scale_val);
    
    /* compute A * phi + b_prime */
    vecScale(_b_prime, -1);
    matMultAdd(_A, _phi_old, _b_prime, _b);
    vecScale(_b_prime, -1);
    sumold = vecSum(_b);

    /* compute keff */
    _k_eff = float(sumnew)/float(sumold);
    log_printf(INFO, "CMFD iter: %i, keff: %.10f, snew: %f, sold: %f, scale_val: %f", iter, _k_eff, sumnew, sumold, scale_val);
    
    /* recompute and normalize initial source */
    matMult(_M, _phi_old, _sold);
    sumold = vecSum(_sold);
    scale_val = (_cells_x*_cells_y * _num_groups) / sumold;
    vecScale(_sold, scale_val);
    sumold = _cells_x*_cells_y * _num_groups;
    
    /* perform power iterations to converge the flux */
    for (iter = 0; iter < 20000; iter++){
      
      /* solve A * b = phi_new */
      vecWAXPY(_b, _k_eff, _b_prime, _sold);
      linearSolve(_A, _phi_new, _b, conv, _omega);
      
      /* computed the new source */
      matMult(_M, _phi_new, _snew);
      sumnew = vecSum(_snew);
      
      /* compute new keff */
      _k_eff = sumnew / sumold;

      /* compute the L2 norm of source error */
      norm = 0.0;
      for (int i = 0; i < _cells_x*_cells_y*_num_groups; i++)
	norm += pow(_snew[i]/_k_eff - _sold[i], 2);
      
      norm = pow(norm, 0.5);
      norm = norm / (_cells_x*_cells_y*_num_groups);
      scale_val = (_cells_x*_cells_y*_num_groups) / sumnew;
      vecScale(_snew, scale_val);
      vecCopy(_snew, _sold);

      log_printf(INFO, "CMFD iter: %i, keff: %f, error: %f", iter + 1, _k_eff, norm);

      /* check for convergence */
      if (norm < _conv_criteria)
	break;
    }
  }
  
  /* rescale flux and pass to meshCells */
  rescaleFlux();
  
  /* give the petsc flux array to the mesh cell flux array */
  setMeshCellFlux();
  
  updateMOCFlux();
  
  return _k_eff;
}


void CmfdTransient::vecWAXPY(double* vec_w, double a, double* vec_x, double* vec_y){

  vecZero(vec_w);

  for (int i = 0; i < _cells_x*_cells_y*_num_groups; i++)
    vec_w[i] = a * vec_x[i] + vec_y[i];

}

void CmfdTransient::matMultAdd(double* mat, double* vec_1, double* vec_2, double* vec_3){

  vecZero(vec_3);

  for (int i = 0; i < _cells_x*_cells_y; i++){
    for (int g = 0; g < _num_groups; g++){
      for (int e = 0; e < _num_groups; e++){
	vec_3[i*_num_groups+g] += mat[(i*_num_groups+g)*_num_groups+e] * vec_1[i*_num_groups+e];
      }
    }
  }

  for (int i = 0; i < _cells_x*_cells_y*_num_groups; i++)
    vec_3[i] += vec_2[i];

}

/* rescale flux */
void CmfdTransient::rescaleFlux(){

  /* initialize variables */
  double sumnew, scale_val;

  /* rescale the new and old flux to have an avg source of 1.0 */
  sumnew = vecSum(_phi_new);
  scale_val = _cells_x*_cells_y*_num_groups / sumnew;
  vecScale(_phi_new, scale_val);
}


/* Fill in the values in the A matrix, M matrix, and phi_old vector
 * @param A matrix
 * @param M matrix
 * @param old flux vector
 * @param solve methed - either DIFFUSION or CMFD
 * @return petsc error indicator
 */
void CmfdTransient::constructMatrices(){

  log_printf(INFO,"Constructing cmfd transient matrices...");
  
  /* initialized variables */
  int row, col;
  double value, phi, b_prime = 0;
  int cell;
  
  Material **materials = _mesh->getMaterials();
  double* old_flux = _mesh->getFluxes(PREVIOUS);
  double* heights = _mesh->getLengthsY();
  double* widths = _mesh->getLengthsX();
  
  vecZero(_phi_old);
  vecZero(_b_prime);
  vecZero(_b);
  matZero(_M, _num_groups);
  matZero(_A, _num_groups*4);

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
	_phi_old[row] = value;
	
	/* absorption term */
	value = materials[cell]->getSigmaA()[e] * _mesh->getVolumes()[cell];
	_A[row*(_num_groups+4)+e+2] += value;

	/* out (1st) and in (2nd) scattering */
	if (_flux_method == PRIMAL){
	  for (int g = 0; g < _num_groups; g++){
	    if (e != g){
	      value = materials[cell]->getSigmaS()[g*_num_groups + e] * _mesh->getVolumes()[cell]; 
	      _A[row*(_num_groups+4)+e+2] += value;	      

	      if (_multigroup_PKE)
		value = - materials[cell]->getSigmaS()[e*_num_groups + g] * _mesh->getVolumes()[cell] * _amplitude[g] / _amplitude[e];
	      else
		value = - materials[cell]->getSigmaS()[e*_num_groups + g] * _mesh->getVolumes()[cell] * _amplitude[0] / _amplitude[0];  

	      _A[row*(_num_groups+4)+g+2] += value;
	    }
	  }
	}
	else{
	  for (int g = 0; g < _num_groups; g++){
	    if (e != g){

	      col = cell*_num_groups+e;
	      value = materials[cell]->getSigmaS()[g*_num_groups + e] * _mesh->getVolumes()[cell];
	      _A[col*(_num_groups+4)+e+2] += value;

	      col = cell*_num_groups+g;
	      value = - materials[cell]->getSigmaS()[e*_num_groups + g] * _mesh->getVolumes()[cell];
	      _A[col*(_num_groups+4)+e+2] += value;
	    }
	  }
	}

	/* add frequency term to diagonal of A */
	if (_transient_method == OMEGA_MODE){
	  if (_multigroup_PKE)
	    value = _frequency[e] / (_velocity[e] * _amplitude[e]) * _mesh->getVolumes()[cell];
	  else
	    value = _frequency[0] / (_velocity[e] * _amplitude[0]) * _mesh->getVolumes()[cell];

	  _A[row*(_num_groups+4)+e+2] += value;
	}
	else if (_transient_method == IQS && !_initial_solve){
	  if (_multigroup_PKE)
	    value = 1.0 / _velocity[e] * (_frequency[e] / _amplitude[e] + 1.0 / _time_stepper->getDtOuter()) * _mesh->getVolumes()[cell];
	  else
	    value = 1.0 / _velocity[e] * (_frequency[0] / _amplitude[0] + 1.0 / _time_stepper->getDtOuter()) * _mesh->getVolumes()[cell];
	  
	  _A[row*(_num_groups+4)+e+2] += value;
	}
	
	/* add previous shape to new vector */
	if (_transient_method == IQS && !_initial_solve){
	  value = 1.0 / (_velocity[e] * _time_stepper->getDtOuter()) * _mesh->getVolumes()[cell] * old_flux[cell*_num_groups + e];
	  _b_prime[row] = value;
	}
		
	/* RIGHT SURFACE */
       
	/* set transport term on diagonal */
	value = (materials[cell]->getDifHat()[2*_num_groups + e] 
			    - materials[cell]->getDifTilde()[2*_num_groups + e]) 
	                    * heights[cell / _cells_x];

	_A[row*(_num_groups+4)+e+2] += value;
	
	/* set transport term on off diagonal */
	if (x != _cells_x - 1){
	  value = - (materials[cell]->getDifHat()[2*_num_groups + e] 
					+ materials[cell]->getDifTilde()[2*_num_groups + e]) 
	                                * heights[cell / _cells_x];

	  _A[row*(_num_groups+4)+_num_groups+2] += value; 
	}

	/* LEFT SURFACE */
	
	/* set transport term on diagonal */
	value = (materials[cell]->getDifHat()[0*_num_groups + e] 
			    + materials[cell]->getDifTilde()[0*_num_groups + e]) 
	                    * heights[cell / _cells_x];
	
	_A[row*(_num_groups+4)+e+2] += value;
	
	/* set transport term on off diagonal */
	if (x != 0){
	    value = - (materials[cell]->getDifHat()[0*_num_groups + e] 
			    - materials[cell]->getDifTilde()[0*_num_groups + e]) 
	                    * heights[cell / _cells_x];

	  _A[row*(_num_groups+4)] += value; 
	}
	
	/* BOTTOM SURFACE */
	
	/* set transport term on diagonal */
	value = (materials[cell]->getDifHat()[1*_num_groups + e] 
			    - materials[cell]->getDifTilde()[1*_num_groups + e]) 
	                    * widths[cell % _cells_x];

	_A[row*(_num_groups+4)+e+2] += value;
       
	/* set transport term on off diagonal */
	if (y != _cells_y - 1){
	  value = - (materials[cell]->getDifHat()[1*_num_groups + e] 
			  + materials[cell]->getDifTilde()[1*_num_groups + e]) 
	                  * widths[cell % _cells_x];

	  _A[row*(_num_groups+4)+1] += value; 
	}
	
	/* TOP SURFACE */
	
	/* set transport term on diagonal */
        value = (materials[cell]->getDifHat()[3*_num_groups + e] 
			    + materials[cell]->getDifTilde()[3*_num_groups + e]) 
	                    * widths[cell % _cells_x];

	_A[row*(_num_groups+4)+e+2] += value;
	
	/* set transport term on off diagonal */
	if (y != 0){
	  value = - (materials[cell]->getDifHat()[3*_num_groups + e] 
					- materials[cell]->getDifTilde()[3*_num_groups + e]) 
	                                * widths[cell % _cells_x];

	  _A[row*(_num_groups+4)+_num_groups+3] += value; 
	}
			
	/* add fission terms to M */
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
	    _M[row*_num_groups+g] += value; 
	  else
	    _M[col*(_num_groups)+e] += value; 
	}
      }
    }
  }
  
  log_printf(INFO,"Done constructing matrices...");
  
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


void CmfdTransient::setInitialSolve(bool initial){
  _initial_solve = initial;
}
