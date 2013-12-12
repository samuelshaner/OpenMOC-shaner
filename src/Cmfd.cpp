#include "Cmfd.h"

/**
 * @brief Constructor initializes boundaries and variables that describe
 *          the cmfd solver object.
 * @details The construcor initializes the many variables that describe
 *          the cmfd mesh, the solve method, and flux type. If solve
 *          method is diffusion theory, the fsr volumes, materials,
 *          and fluxes are initialized.
 * @param geometry pointer to the geometry
 * @param criteria convergence criteria on keff
 */
Cmfd::Cmfd(Geometry* geometry, double criteria) {

    /* Create objects */
    _geometry = geometry;
    _mesh = geometry->getMesh();
    _timer = new Timer();

    /* integer values */
    _cx = _mesh->getCellsX();
    _cy = _mesh->getCellsY();
    _ng = _mesh->getNumGroups();
    _nc = _ng*_cx*_cy;
    
    /* float values */
    _omega = 1.0;
    _conv_criteria = criteria;
    
    /* Boolean and Enum flags to toggle features */
    _solve_method = _mesh->getSolveType();
    
    /* Create Cmfd matrix objects */
    _M        = new double[_nc*_ng];
    _A        = new double[_nc*(4+_ng)];
    _phi_temp = new double[_nc];  
    _snew     = new double[_nc];  
    _sold     = new double[_nc];
    _phi_old  = NULL;  
    _phi_new  = NULL;  
    _b        = new double[_nc];
    _b_prime  = new double[_nc];  

    _AM       = new double[_nc*(4+_ng)];
    _y        = new double[_nc];
    
    /* If solving diffusion problem, create arrays for FSR parameters */
    if (_solve_method == DIFFUSION)
	_mesh->initializeSurfaceCurrents();
}


/**
 * @Brief Destructor deletes arrays of A and M row insertion arrays
 */
Cmfd::~Cmfd() {

  delete [] _phi_temp;
  delete [] _A;
  delete [] _M;
  delete [] _snew;
  delete [] _sold;
  delete [] _y;
  delete [] _AM;
  delete [] _b;
  delete [] _b_prime;

}


/*
 * @brief CMFD solver that solves the diffusion problem
 * @return k-effective
 */
double Cmfd::computeKeff(){

    log_printf(INFO, "Running CMFD diffusion solver...");
    
    /* if solving diffusion problem, initialize timer */
    if (_solve_method == DIFFUSION)
	_timer->startTimer();
    
    /* initialize variables */
    double sumnew = 0.0;
    double sumold = 0.0;
    double norm = 1e10;
    double scale_val;
    double conv = 1e-8;
    int iter = 0;
    int method = 1;
    double balance;

    if (_solve_method == MOC)
	_mesh->computeXS();
    else
	_timer->startTimer();

    _mesh->computeDs();

    _phi_old = _mesh->getFluxes(FSR_OLD);
    _mesh->copyFlux(FSR_OLD, CURRENT);
    _phi_new = _mesh->getFluxes(CURRENT);
    
    constructMatrices();

    if (_mesh->getInitialState() == false){
      matSubtract(_AM, _A, 1.0, _M, _cx, _cy, _ng);
      linearSolveRB(_AM, _phi_new, _b, _phi_temp, _conv_criteria, _omega, _cx, _cy, _ng, 10000);

    }
    else{
	if (method == 1){
	    
	    /* get initial source */
	    matMultM(_M, _phi_old, _snew, _cx*_cy, _ng);
	    
	    /* scale the initial source */
	    sumnew = vecSum(_snew, _nc);
	    scale_val = _nc / sumnew;
	    vecScale(_snew, scale_val, _nc);
	    vecCopy(_phi_old, _phi_new, _nc);		
	    sumold = _nc;
	    	    
	    vecCopy(_snew, _sold, _nc);
	    
	    while (norm > _conv_criteria){
		
		/* Solver phi = A^-1 * old_source */
		vecWAXPY(_b_prime, 1.0, _snew, _b, _nc);
		linearSolveRB(_A, _phi_new, _b_prime, _phi_temp, conv, _omega, _cx, _cy, _ng, 1000);
		
		/* compute new source */
		matMultM(_M, _phi_new, _snew, _cx*_cy, _ng);
		
		/* compute keff */
		sumnew = vecSum(_snew, _nc);
		_k_eff = sumnew / sumold;
		vecScale(_sold, _k_eff, _nc);
				
		/* compute l2 norm */
		norm = 0.0;
		for (int i = 0; i < _nc; i++)
		    if (_snew[i] != 0.0)
			norm += pow((_snew[i] - _sold[i])/_snew[i], 2);
	    
		norm = pow(norm, 0.5);
		norm = norm / _nc;
		
		/* scale the new source and pass to old source */
		scale_val = _nc / sumnew;
		vecScale(_snew, scale_val, _nc);
		
		vecCopy(_snew, _sold, _nc);
		
		iter++;

		log_printf(INFO, "iter: %i, k_eff: %f, error: %f", iter, _k_eff, norm);
		
		if (iter == 10000)
		    break;
	    }	
	    
	    /* converge the flux with 1000 GS iterations */
	    //vecWAXPY(_b_prime, 1.0, _snew, _b, _nc);
	    //linearSolve(_A, _phi_new, _b_prime, _phi_temp, -1.0, _omega, _cx, _cy, _ng, 1000);
	    
	}
	else{

	    log_printf(INFO, "solving eigenvalue problem with wielandt shift method");
	    
	    int max_iter = 100;
	    double shift, offset;
	    
	    /* copy old flux to new */
	    vecCopy(_phi_old, _phi_new, _nc);
	    
	    /* compute the initial rayleigh coefficient */
	    _k_eff = rayleighQuotient(_A, _M, _phi_new, _snew, _phi_temp, _cx, _cy, _ng);
	    
	    /* compute the initial old source */
	    matMultM(_M, _phi_new, _sold, _cx*_cy, _ng);
	    sumold = vecSum(_sold, _nc);
	    scale_val = _nc / sumold;
	    vecScale(_sold, scale_val, _nc);
	    vecScale(_phi_new, scale_val, _nc);
	    sumold = _nc;
	    
	    vecCopy(_phi_new, _y, _nc);
	    
	    log_printf(INFO, "initial keff: %f", _k_eff);
	    
	    shift = 1.5;

	    /* loosely converge eigenvaue with shift of 1.5 */
	    while (norm > 1e-5){
		
		/* reconstruct _AM */
		matSubtract(_AM, _A, 1.0/shift, _M, _cx, _cy, _ng);
		
		/* solve inverse system */
		linearSolveRB(_AM, _y, _snew, _phi_temp, conv, _omega, _cx, _cy, _ng, max_iter);
		
		/* compute new _phi_new  */
		vecCopy(_y, _phi_new, _nc);
		vecScale(_phi_new, vecMax(_y, _nc), _nc);
		matMultM(_M, _phi_new, _snew, _cx*_cy, _ng);
		sumnew = vecSum(_snew, _nc);
		scale_val = _nc / sumnew;
		vecScale(_snew, scale_val, _nc);
		vecScale(_phi_new, scale_val, _nc);
		
		/* compute new eigenvalue */
		_k_eff = rayleighQuotient(_A, _M, _phi_new, _snew, _phi_temp, _cx, _cy, _ng);
		
		/* compute norm */
		norm = 0.0;
		for (int i = 0; i < _nc; i++)
		    norm += pow(_snew[i] - _sold[i], 2);
		
		norm = pow(norm, 0.5);
		norm = norm / _nc;
		
		/* copy new source to old */
		vecCopy(_snew, _sold, _nc);
		
		iter++;
		
		log_printf(INFO, "iter: %i, _k_eff: %f, norm: %f", iter, _k_eff, norm); 
	    }

	    offset = 0.05;
	    log_printf(INFO, "setting offset to %f", offset);

	    /* finely converge eigenvalue with constant offset */
	    while (norm > _conv_criteria){
		
		/* reconstruct _AM */
	    	matSubtract(_AM, _A, 1.0/(_k_eff + offset), _M, _cx, _cy, _ng);
		
		/* solve inverse system */
	    	linearSolveRB(_AM, _y, _snew, _phi_temp, conv, _omega, _cx, _cy, _ng, max_iter);
		
		/* compute new _phi_new  */
	    	vecCopy(_y, _phi_new, _nc);
		vecScale(_phi_new, vecMax(_y, _nc), _nc);
		matMultM(_M, _phi_new, _snew, _cx*_cy, _ng);
		sumnew = vecSum(_snew, _nc);
		scale_val = _nc / sumnew;
		vecScale(_snew, scale_val, _nc);
		vecScale(_phi_new, scale_val, _nc);
		
		/* compute new eigenvalue */
		_k_eff = rayleighQuotient(_A, _M, _phi_new, _snew, _phi_temp, _cx, _cy, _ng);
		
		/* compute norm */
		norm = 0.0;
		for (int i = 0; i < _nc; i++)
		    norm += pow(_snew[i] - _sold[i], 2);
		
		norm = pow(norm, 0.5);
		norm = norm / _nc;
		
		/* copy new source to old */
		vecCopy(_snew, _sold, _nc);
		
		iter++;
		
		log_printf(INFO, "iter: %i, _k_eff: %f, norm: %f", iter, _k_eff, norm); 
	    }
	}
    }
	
    
    /* rescale the old and new flux */
    if (_mesh->getInitialState() == true && _solve_method == MOC)
	rescaleFlux();
    
    /* update the MOC flux */
    if (_solve_method == MOC)
	_mesh->updateMOCFlux();  
    
    /* If solving diffusion problem, print timing results */
    if (_solve_method == DIFFUSION){
	std::string msg_string;
	log_printf(TITLE, "TIMING REPORT");
	_timer->stopTimer();
	_timer->recordSplit("Total time to solve diffusion eigenvalue problem");
	
	double tot_time = _timer->getSplit("Total time to solve diffusion eigenvalue problem");
	msg_string = "Total time to solve diffusion eigenvalue problem";
	msg_string.resize(53, '.');
	log_printf(RESULT, "%s%1.4E sec", msg_string.c_str(), tot_time);
    }

    return _k_eff;
}


/**
 * @brief rescale the initial and converged flux arrays
 * @return petsc_err petsc error flag
 */
void Cmfd::rescaleFlux(){

    double sumnew, sumold, scale_val;
    
    matMultM(_M, _phi_new, _snew, _cx*_cy, _ng);
    sumnew = vecSum(_snew, _nc);
    scale_val = _nc / sumnew;
    vecScale(_phi_new, scale_val, _nc);
    matMultM(_M, _phi_old, _sold, _cx*_cy, _ng);
    sumold = vecSum(_sold, _nc);
    scale_val = _nc / sumold;
    vecScale(_phi_old, scale_val, _nc);
    
}


/* Fill in the values in the A matrix, M matrix, and phi_old vector
 * @return petsc_err petsc error flag
 */
void Cmfd::constructMatrices(){

    log_printf(INFO,"Constructing matrices...");
    
    double value;
    int cell;
    int row;
    int offset = int(FORWARD)*4*_ng;
    int e, g, dg;
    
    /* get arrays */
    Material** materials = _mesh->getMaterials();
    Material* material;
    double* heights = _mesh->getLengthsY();
    double* widths = _mesh->getLengthsX();
    double* volumes = _mesh->getVolumes();
    double dt;
    double* velocity;
    double* frequency_new;
    double* flux_old;
    TimeStepper* ts;

    if (_mesh->getInitialState() == false){
	dt = _mesh->getDtMOC();
	velocity = _mesh->getVelocity();
	frequency_new = _mesh->getFrequencies(CURRENT);
	flux_old = _mesh->getFluxes(PREVIOUS_CONV);
	ts = _mesh->getTimeStepper();
    }

    vecZero(_b, _nc);
    vecZero(_M, _ng*_nc);
    vecZero(_A, (_ng+4)*_nc);
    vecZero(_AM, (_ng+4)*_nc);
    
    /* loop over cells */
    #pragma(cell, row, value, material)
    for (int y = 0; y < _cy; y++){
	for (int x = 0; x < _cx; x++){
	    
	    cell = y*_cx + x;
	    material = materials[cell];

	    /* loop over groups */
	    for (e = 0; e < _ng; e++){
		
		row = cell*_ng + e;

		if (_mesh->getInitialState() == false){
		    
		    /* method source */
		    if (_mesh->getTransientType() == THETA){

			_b[row] = flux_old[row] / (velocity[e] * dt) * volumes[cell];

			value = volumes[cell] / (velocity[e] * dt); 
			_A[row*(_ng+4)+e+2] += value;
		    }
		    else if (_mesh->getTransientType() == MAF){

		        _b[row] = _mesh->getFluxes(PREVIOUS)[row] 
			  / (velocity[e] * dt) * volumes[cell]; 

			value = volumes[cell] / velocity[e] * (1.0 / dt + frequency_new[row]); 
			_A[row*(_ng+4)+e+2] += value;
		    }
		    else if (_mesh->getTransientType() == ADIABATIC){

			_b[row] = flux_old[row] / (velocity[e] * dt) *
			  exp(frequency_new[row] * dt) 
			  * volumes[cell]; 
	
			value = volumes[cell] / (velocity[e] * dt) 
			    * (1.0 + log(_mesh->getFluxes(PREVIOUS)[row] / flux_old[row])); 
			_A[row*(_ng+4)+e+2] += value;

		    }

		    /* delayed source */
		    if (material->getType() == FUNCTIONAL){
			for (dg = 0; dg < _mesh->getNumDelayGroups(); dg++){
			    _b[row] += material->getChi()[e] * volumes[cell] * 
				_mesh->getLambda()[dg] * static_cast<FunctionalMaterial*>(material)->getPrecConc(CURRENT, dg);
			}
		    }	
		}	
				
		/* absorption term */
		value = material->getSigmaA()[e] * volumes[cell];
		_A[row*(_ng+4)+e+2] += value;

		/* out (1st) and in (2nd) scattering */
		for (g = 0; g < _ng; g++){
		    if (e != g){
			value = material->getSigmaS()[g*_ng+e] * volumes[cell]; 
			_A[row*(_ng+4)+e+2] += value;
			value = - material->getSigmaS()[e*_ng + g] * volumes[cell];
			_A[row*(_ng+4)+g+2] += value;
		    }
		}
				
		/* RIGHT SURFACE */
		
		/* set transport term on diagonal */
		
		value = (material->getDifHat()[2*_ng + e] 
			 - material->getDifTilde()[offset + 2*_ng + e]) 
		  * heights[cell / _cx];
		
		_A[row*(_ng+4)+e+2] += value; 
		
		
		/* set transport term on off diagonal */
		if (x != _cx - 1){
		    value = - (material->getDifHat()[2*_ng + e] 
			       + material->getDifTilde()[offset + 2*_ng + e]) 
		      * heights[cell / _cx];
		    
		    _A[row*(_ng+4)+_ng+2] += value; 
		}
		
		/* LEFT SURFACE */
		
		/* set transport term on diagonal */
		value = (material->getDifHat()[0*_ng + e] 
			 + material->getDifTilde()[offset + 0*_ng + e]) 
		    * heights[cell / _cx];
		
		_A[row*(_ng+4)+e+2] += value; 
		
		/* set transport term on off diagonal */
		if (x != 0){
		    value = - (material->getDifHat()[0*_ng + e] 
			       - material->getDifTilde()[offset + 0*_ng + e]) 
		      * heights[cell / _cx];
		    
		    _A[row*(_ng+4)] += value; 
		}
		
		/* BOTTOM SURFACE */
		
		/* set transport term on diagonal */
		value = (material->getDifHat()[1*_ng + e] 
			 - material->getDifTilde()[offset + 1*_ng + e]) 
		  * widths[cell % _cx];
		
		_A[row*(_ng+4)+e+2] += value;        
		
		/* set transport term on off diagonal */
		if (y != _cy - 1){
		    value = - (material->getDifHat()[1*_ng + e] 
			       + material->getDifTilde()[offset + 1*_ng + e]) 
		      * widths[cell % _cx];
		    
		    _A[row*(_ng+4)+1] += value; 
		}
		
		/* TOP SURFACE */
		
		/* set transport term on diagonal */
		value = (material->getDifHat()[3*_ng + e] 
			 + material->getDifTilde()[offset + 3*_ng + e]) 
		  * widths[cell % _cx];
		
		_A[row*(_ng+4)+e+2] += value; 
		
		/* set transport term on off diagonal */
		if (y != 0){
		    value = - (material->getDifHat()[3*_ng + e] 
			       - material->getDifTilde()[offset + 3*_ng + e]) 
		      * widths[cell % _cx];
		    
		    _A[row*(_ng+4)+_ng+3] += value; 
		}
		
		/* fission source */
		for (g = 0; g < _ng; g++){	    
		    if (_mesh->getInitialState() == true){
			value = material->getChi()[e] * material->getNuSigmaF()[g] *  volumes[cell];
		    }
		    else{
			value = (1.0 - _mesh->getBetaSum()) / _mesh->getKeff0() * 
			    material->getChi()[e] * material->getNuSigmaF()[g] * 
			    volumes[cell];	    
		    }
		    
		    _M[row*_ng+g] += value; 
		}
	    }
	}
    }

    log_printf(INFO,"Done constructing matrices...");  
    
}


/**
 * @brief get pointer to loss matrix, A
 * @return _A pointer to loss matrix, A
 */
double* Cmfd::getA(){
    return _A;
}


/**
 * @brief get pointer to source matrix, M
 * @return _M pointer to source matrix, M
 */
double* Cmfd::getM(){
    return _M;
}

/**
 * @brief get k_eff
 * @return _k_eff k_eff
 */
double Cmfd::getKeff(){
    return _k_eff;
}


/**
 * @brief Get pointer to the mesh 
 * @return _mesh pointer to mesh
 */
Mesh* Cmfd::getMesh(){
  return _mesh;
}


void Cmfd::setOmega(double omega){
  _omega = omega;
}


void Cmfd::checkNeutronBalance(){

  //    if (_solve_method == MOC)
  //	_mesh->computeXS();

  //_mesh->computeDs();

    if (_solve_method == MOC)
        _phi_new = _mesh->getFluxes(FSR_OLD);
    else
        _phi_new = _mesh->getFluxes(CURRENT);
    
    /* construct matrices */
    constructMatrices();
    
    /* get right hand side */
    matMultM(_M, _phi_new, _snew, _cx*_cy, _ng);
    vecScale(_snew, 1.0/_mesh->getKeff0(), _nc);
    
    matMultA(_A, _phi_new, _sold, _cx, _cy, _ng);
    
    double sumleft = vecSum(_sold, _nc);
    double sumright = vecSum(_snew, _nc);
    
    log_printf(NORMAL, "cmfd balance: %1.3E", sumleft - sumright);
    
}


void Cmfd::setNumThreads(int num_threads) {

    if (num_threads <= 0)
        log_printf(ERROR, "Unable to set the number of threads for the Solver "
		   "to %d since it is less than or equal to 0", num_threads);

    /* Set the number of threads for OpenMP */
    omp_set_num_threads(num_threads);
}
