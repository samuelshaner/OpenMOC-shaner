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
    _cells_x = _mesh->getCellsX();
    _cells_y = _mesh->getCellsY();
    _num_groups = _mesh->getNumGroups();
    _nc = _num_groups*_cells_x*_cells_y;
    
    /* float values */
    _omega = 1.0;
    _conv_criteria = criteria;
    
    /* Boolean and Enum flags to toggle features */
    _solve_method = _mesh->getSolveType();
    
    /* Create Cmfd matrix objects */
    _M        = new double[_nc*_num_groups];
    _A        = new double[_nc*(4+_num_groups)];
    _phi_temp = new double[_nc];  
    _snew     = new double[_nc];  
    _sold     = new double[_nc];
    _phi_old  = NULL;  
    _phi_new  = NULL;  
    
    /* If solving diffusion problem, create arrays for FSR parameters */
    if (_solve_method == DIFFUSION)
	_mesh->initializeSurfaceCurrents();
}


/**
 * @brief Destructor deletes arrays of A and M row insertion arrays
 */
Cmfd::~Cmfd() {

  delete [] _phi_temp;
  delete [] _A;
  delete [] _M;
  delete [] _snew;
  delete [] _sold;

}


/*
 * @brief CMFD solver that solves the diffusion problem
 * @return k-effective
 */
double Cmfd::computeKeff(){

    log_printf(NORMAL, "Running CMFD diffusion solver...");
    
    /* if solving diffusion problem, initialize timer */
    if (_solve_method == DIFFUSION)
	_timer->startTimer();
    
    /* initialize variables */
    double sumnew = 0.0;
    double sumold = 0.0;
    double norm = 0.0;
    double scale_val;
    double conv = 1e-2;

    if (_solve_method == MOC)
	_mesh->computeXS();
    
    _mesh->computeDs();

    _phi_old = _mesh->getFluxes(PREVIOUS);
    _phi_new = _mesh->getFluxes(CURRENT);
    
    constructMatrices();
    
    /* get initial source */
    matMultM(_M, _phi_old, _sold, _cells_x*_cells_y, _num_groups);
    sumold = vecSum(_sold, _nc);
    scale_val = _nc / sumold;
    vecScale(_sold, scale_val, _nc);
    
    sumold = _nc;
    
    for (int iter = 0; iter < 20000; iter++){
	
	/* solver phi = A^-1 * old_source */
	linearSolve(_A, _phi_new, _sold, _phi_temp, conv, _omega, _cells_x, _cells_y, _num_groups);
	
	/* compute new source */
	matMultM(_M, _phi_new, _snew, _cells_x*_cells_y, _num_groups);
	sumnew = vecSum(_snew, _nc);
	
	/* compute keff */
	_k_eff = sumnew / sumold;
	
	vecScale(_sold, _k_eff, _nc);
	
	/* compute l2 norm */
	norm = 0.0;
	for (int i = 0; i < _nc; i++)
	    norm += pow(_snew[i] - _sold[i], 2);
	
	norm = pow(norm, 0.5);
	norm = norm / _nc;

	/* scale the new source and pass to old source */
	scale_val = _nc / sumnew;
	vecScale(_snew, scale_val, _nc);
	vecCopy(_snew, _sold, _nc);
	
	log_printf(INFO, "GS POWER iter: %i, keff: %f, error: %f", iter, _k_eff, norm);
	
	if (norm < _conv_criteria)
	    break;
    }

    /* rescale the old and new flux */
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
    
    matMultM(_M, _phi_new, _snew, _cells_x*_cells_y, _num_groups);
    sumnew = vecSum(_snew, _nc);
    scale_val = _nc / sumnew;
    vecScale(_phi_new, scale_val, _nc);
    matMultM(_M, _phi_old, _sold, _cells_x*_cells_y, _num_groups);
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
    
    /* get arrays */
    Material** materials = _mesh->getMaterials();
    double* heights = _mesh->getLengthsY();
    double* widths = _mesh->getLengthsX();
    
    matZero(_M, _num_groups, _nc);
    matZero(_A, _num_groups+4, _nc);
    
    /* loop over cells */
    for (int y = 0; y < _cells_y; y++){
	for (int x = 0; x < _cells_x; x++){
	    
	    cell = y*_cells_x + x;
	    
	    /* loop over groups */
	    for (int e = 0; e < _num_groups; e++){
		
		row = cell*_num_groups + e;
		
		/* absorption term */
		value = materials[cell]->getSigmaA()[e] * _mesh->getVolumes()[cell];
		_A[row*(_num_groups+4)+e+2] += value;
		
		/* out (1st) and in (2nd) scattering */
		for (int g = 0; g < _num_groups; g++){
		    if (e != g){
			value = materials[cell]->getSigmaS()[g*_num_groups+e] * _mesh->getVolumes()[cell]; 
			_A[row*(_num_groups+4)+e+2] += value;
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
		
		/* source term */
		for (int g = 0; g < _num_groups; g++){	    
		    value = materials[cell]->getChi()[e] * materials[cell]->getNuSigmaF()[g] * _mesh->getVolumes()[cell];	    
		    
		    _M[row*_num_groups+g] += value; 
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
