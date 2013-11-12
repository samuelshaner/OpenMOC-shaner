/*
 * Tcmfd.cpp
 *
 *  Created on: September 8, 2013
 *  Author: Sam Shaner
 *  MIT, Course 22
 *  shaner@mit.edu
 */

#include "TCmfd.h"

/**
 * Tcmfd constructor
 * @param geom pointer to the geometry
 * @param criteria convergence criterial on keff 
 */
Tcmfd::Tcmfd(Geometry* geometry, double criteria){
    
    /* Boolean and Enum flags to toggle features */
    _transient_method = ADIABATIC;
    _geom  = geometry;
    _mesh  = geometry->getMesh();
    _timer = new Timer();

    /* integert values */
    _cells_x    = _mesh->getCellsX();
    _cells_y    = _mesh->getCellsY();
    _num_groups = _mesh->getNumGroups();
    _nc         = _num_groups*_cells_x*_cells_y;

    /* float values */
    _omega         = 1.0;
    _conv_criteria = criteria;
    _k_eff_0       = 1.0;
    _beta_sum      = 0.0;
    _dt_cmfd       = 1e-4;

    /* Boolean and Enum flags to toggle features */
    _solve_method = _mesh->getSolveType();
    
    /* Create matrix and vector objects */
    _M         = new double[_nc*_num_groups];
    _A         = new double[_nc*(4+_num_groups)];
    _phi_temp  = new double[_nc];  
    _snew      = new double[_nc];  
    _sold      = new double[_nc];
    _b         = new double[_nc];
    _b_prime   = new double[_nc];
    _phi_new   = NULL;
    _phi_old   = NULL;

    /* transient parameters */
    _lambda        = NULL;
    _beta          = NULL;
    _velocity      = NULL;
    _initial_state = true;

    /* If solving diffusion problem, create arrays for FSR parameters */
    if (_solve_method == DIFFUSION)
	_mesh->initializeSurfaceCurrents();
    
}

/**
 * cmfd Destructor clears all memory
 */
Tcmfd::~Tcmfd() {

    delete [] _b;
    delete [] _b_prime;
    delete [] _M;
    delete [] _A;
    delete [] _phi_temp;
    delete [] _snew;
    delete [] _sold;

}



/*
 * CMFD solver that solves the diffusion problem
 * @param solve methed - either diffusion or cmfd (acceleration)
 * @param iteration number of in MOC solver - used for plotting
 * @return k-effective
 */
void Tcmfd::solveTCMFD(){

    log_printf(INFO, "Running cmfd diffusion diffusion solver...");
    
    /* initialize variables */
    int iter = 0;
    double norm = 0.0;
    double conv = 1e-3;
    
    if (_solve_method == MOC)
	_mesh->computeXS();
    
    _mesh->computeDs();
    
    _phi_old = _mesh->getFluxes(PREVIOUS);
    _phi_new = _mesh->getFluxes(CURRENT);

    /* construct matrices */
    constructMatrices();
    
    /* get initial source and find initial k_eff */
    matMultM(_M, _phi_new, _snew, _cells_x*_cells_y, _num_groups);

    /* perform power iterations to converge the flux */
    for (iter = 0; iter < 20000; iter++){
	
	/* solve A * b = phi_new */
	vecWAXPY(_b_prime, 1.0, _snew, _b, _nc);
	linearSolve(_A, _phi_new, _b_prime, _phi_temp, conv, _omega, _cells_x, _cells_y, _num_groups);
	
	/* computed the new source */
	matMultM(_M, _phi_new, _snew, _cells_x*_cells_y, _num_groups);
	
	/* compute the L2 norm of source error */
	norm = 0.0;
	for (int i = 0; i < _nc; i++)
	    norm += pow((_snew[i] - _sold[i])/(_snew[i]+1e-15), 2);

	norm = pow(norm, 0.5);
	norm = norm / _nc;

	/* copy new source to old */
	vecCopy(_snew, _sold, _nc);
	
	/* check for convergence */
	if (norm < _conv_criteria)
	    break;    
    }

    log_printf(INFO, "CMFD iter: %i", iter);
}


/* Fill in the values in the A matrix, M matrix, and phi_old vector
 * @param A matrix
 * @param M matrix
 * @param old flux vector
 * @param solve methed - either DIFFUSION or CMFD
 * @return petsc error indicator
 */
void Tcmfd::constructMatrices(){

    log_printf(INFO,"Constructing cmfd transient matrices transient...");
    
    /* initialized variables */
    int row, col;
    double value;
    int cell;
    
    Material **materials = _mesh->getMaterials();
    double* heights = _mesh->getLengthsY();
    double* widths = _mesh->getLengthsX();
    
    vecZero(_b, _nc);
    matZero(_M, _num_groups, _nc);
    matZero(_A, _num_groups*4, _nc);

    /* loop over mesh cells in y direction */
    for (int y = 0; y < _cells_y; y++){
    
	/* loop over mesh cells in x direction */
	for (int x = 0; x < _cells_x; x++){
	    
	    cell = y*_cells_x + x;
	    
	    /* loop over energy groups */
	    for (int e = 0; e < _num_groups; e++){
		
		row = cell*_num_groups+e;
		
		/* old flux and delayed neutron precursors */
		_b[row] = _phi_old[cell*_num_groups + e] / (_velocity[e] * _dt_cmfd) * _mesh->getVolumes()[cell];

		if (materials[cell]->getType() == FUNCTIONAL){
		    for (int dg = 0; dg < _num_delay_groups; dg++){
			_b[row] += materials[cell]->getChi()[e] * _mesh->getVolumes()[cell] * _lambda[dg] * static_cast<FunctionalMaterial*>(materials[cell])->getPrecConc(CURRENT, dg);
		    }
		}		
		
		/* absorption term */
		value = materials[cell]->getSigmaA()[e] * _mesh->getVolumes()[cell];
		_A[row*(_num_groups+4)+e+2] += value;
		
		/* flux derivative term */
		if (!_initial_state){
		    value = _mesh->getVolumes()[cell] / (_velocity[e] * _dt_cmfd);
		    _A[row*(_num_groups+4)+e+2] += value;
		}
		
		/* scattering terms */
		for (int g = 0; g < _num_groups; g++){
		    if (e != g){
			/* out scattering */
			value = materials[cell]->getSigmaS()[g*_num_groups + e] * _mesh->getVolumes()[cell]; 
			_A[row*(_num_groups+4)+e+2] += value;	      
			
			/* in scattering */
			value = - materials[cell]->getSigmaS()[e*_num_groups + g] * _mesh->getVolumes()[cell];			    			    
			_A[row*(_num_groups+4)+g+2] += value;
		    }
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
		    
		    if (_initial_state)
			value = materials[cell]->getChi()[e] * materials[cell]->getNuSigmaF()[g] * _mesh->getVolumes()[cell]/_k_eff_0;
		    else
			value = (1.0 - _beta_sum) * materials[cell]->getChi()[e] * materials[cell]->getNuSigmaF()[g] * _mesh->getVolumes()[cell] / _k_eff_0;
		    			
		    col = cell*_num_groups+g;
		    
		    _M[row*_num_groups+g] += value; 
		}
	    }
	}
    }
    
    log_printf(INFO,"Done constructing matrices...");
}


void Tcmfd::setTransientType(transientType trans_type){
  _transient_method = trans_type;
}


transientType Tcmfd::getTransientType(){
  return _transient_method;
}


void Tcmfd::setKeff0(double keff_0){
  _k_eff_0 = keff_0;
}


void Tcmfd::setTimeStepper(TimeStepper* time_stepper){
    _time_stepper = time_stepper; 
}


double* Tcmfd::getBeta(){
  return _beta;
}


double* Tcmfd::getLambda(){
  return _lambda;
}


double* Tcmfd::getVelocity(){
  return _velocity;
}


void Tcmfd::setBeta(double* beta, int num_delay_groups){

  _beta = new double[num_delay_groups];
  _beta_sum = 0.0;

  for (int g = 0; g < num_delay_groups; g++){
    _beta[g] = beta[g];
    _beta_sum += beta[g];
    log_printf(NORMAL, "g: %i, beta: %f, beta_sum: %f", g, beta[g], _beta_sum);
  }  
}


void Tcmfd::setLambda(double* decay_const, int num_delay_groups){

  _lambda = new double[num_delay_groups];
  
  for (int g = 0; g < num_delay_groups; g++)
    _lambda[g] = decay_const[g];  

}


void Tcmfd::setVelocity(double* velocity, int num_groups){

  _velocity = new double[num_groups];
  
  for (int g = 0; g < num_groups; g++)
    _velocity[g] = velocity[g];
  
}


void Tcmfd::setDtCMFD(double dt){
    _dt_cmfd = dt;
}


void Tcmfd::setInitialState(bool state){
    _initial_state = state;
}


void Tcmfd::setNumDelayGroups(int num_groups){
    _num_delay_groups = num_groups;
}


Mesh* Tcmfd::getMesh(){
    return _mesh;
}


void Tcmfd::setOmega(double omega){
    _omega = omega;
}
