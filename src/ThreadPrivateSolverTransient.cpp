#include "ThreadPrivateSolverTransient.h"


/**
 * @brief Constructor initialize array pointers for tracks and materials.
 * @details The constructor retrieves the number of energy groups and flat
 *          source regions and azimuthal angles from the geometry and track
 *          generator. The constructor initalizes the number of threads to a 
 *          default of 1.
 * @param geometry an optional pointer to the geometry
 * @param track_generator an optional pointer to the trackgenerator
 */
ThreadPrivateSolverTransient::ThreadPrivateSolverTransient(Geometry* geometry, 
					 TrackGenerator* track_generator,
					 Cmfd* cmfd) :
  CPUSolver(geometry, track_generator, cmfd) {

    _thread_flux = NULL;
    _thread_currents = NULL;
    _transient_method = MAF;
}


/**
 * @brief Destructor calls Solver subclass destructor to deletes arrays
 *        for fluxes and sources.
 */
ThreadPrivateSolverTransient::~ThreadPrivateSolverTransient() { 

    if (_thread_flux != NULL) {
        delete [] _thread_flux;
	_thread_flux = NULL;
    }

    if (_thread_currents != NULL) {
        delete [] _thread_currents;
	_thread_currents = NULL;
    }

}


void ThreadPrivateSolverTransient::initialize(){

    /* An initial guess for the eigenvalue */
    _k_eff = 1.0;

    /* Initialize data structures */
    initializePolarQuadrature();
    initializeFluxArrays();
    initializeSourceArrays();
    precomputePrefactors();
    initializeFSRs();

    /* Check that each FSR has at least one segment crossing it */
    checkTrackSpacing();

    /* initialize cmfd */
    initializeCmfd();

    /* Set scalar flux to unity for each region */
    flattenFSRFluxes(1.0);
    flattenFSRSources(1.0);
    zeroTrackFluxes();

}


void ThreadPrivateSolverTransient::resetSegmentMaterials(){

    log_printf(NORMAL, "resetting segment materials");

    /* reset the material pointer for each segment */
    int tid;
    Track* curr_track;
    int num_segments;
    segment* curr_segment;    
    segment* segments;

    /* Loop over azimuthal angle halfspaces */
    for (int i=0; i < 2; i++) {
	
        /* Compute the minimum and maximum track IDs corresponding to 
         * this azimuthal angular halfspace */
        int min = i * (_tot_num_tracks / 2);
	int max = (i + 1) * (_tot_num_tracks / 2);
	
	/* Loop over each thread within this azimuthal angle halfspace */
	#pragma omp parallel for private(tid, curr_track, \
	  num_segments, segments, curr_segment) schedule(guided)
	for (int track_id=min; track_id < max; track_id++) {

	    tid = omp_get_thread_num();
	    
	    /* Initialize local pointers to important data structures */	
	    curr_track = _tracks[track_id];
	    num_segments = curr_track->getNumSegments();
	    segments = curr_track->getSegments();
	    
	    /* Loop over each segment in forward direction */
	    for (int s=0; s < num_segments; s++) {
	        curr_segment = &segments[s];
		curr_segment->_material = _FSR_materials[curr_segment->_region_id];
	    }
        }
    }
}


/**
 * @brief Computes keff by performing a series of transport sweep and 
 *        source updates.
 * @details This is the main method exposed to the user through the Python
 *          interface to run a simulation. The method makes an initial guess
 *          for the scalar and boundary fluxes and peforms transport sweeps
 *          and source updates until convergence.
 * @param max_iterations the maximum number of iterations allowed
 * @return the value of keff computed
 */
FP_PRECISION ThreadPrivateSolverTransient::convergeSource(int max_iterations) {

    log_printf(INFO, "converging transient source");

    _timer->startTimer();

    /* Counter for the number of iterations to converge the source */
    _num_iterations = 0;

    Mesh* mesh = _geometry->getMesh();

    double cmfd_k_eff;

    /* The residual on the source */
    FP_PRECISION residual = 0.0;

    if (_geometry->getMesh()->getInitialState())
      log_printf(NORMAL, "Iter %d: \tk_eff = %1.6f \tCMFD k_eff = %1.6f"
		 "\tres = %1.3E", 0, _k_eff, _k_eff, residual);
    else
      log_printf(NORMAL, "Iter %d: \trho = $%1.4f"
		 "\tres = %1.3E", 0, getReactivity(), residual);


    /* Fission source iteration loop */
    for (int i=0; i < max_iterations; i++) {

	/* normalize the fluxes */
        if (_geometry->getMesh()->getInitialState()){
	    _timer->startTimer();
	    normalizeFluxes();
	    _timer->stopTimer();
	    _timer->recordSplit("MOC normalize fluxes");
	}

	/* compute the new source */
	_timer->startTimer();
	residual = computeFSRSources();
	_timer->stopTimer();
	_timer->recordSplit("MOC compute sources");
	
	/* perform transport sweep */
	_timer->startTimer();
	transportSweep();
	_timer->stopTimer();
	_timer->recordSplit("MOC transport sweep");
	
	/* add source to scalar flux */
	_timer->startTimer();
	addSourceToScalarFlux();
	_timer->stopTimer();
	_timer->recordSplit("MOC add source to flux");
	
	/* update the flux with cmfd */
	if (!(i >= 1 && residual < _source_convergence_thresh)){    
	    _timer->startTimer();
	    cmfd_k_eff = _cmfd->computeKeff();
	    _timer->stopTimer();
	    _timer->recordSplit("CMFD solve");
	}
	
	_timer->startTimer();
	computeKeff();
	_timer->stopTimer();
	_timer->recordSplit("MOC compute keff");

	_num_iterations++;

	if (_geometry->getMesh()->getInitialState())
	  log_printf(NORMAL, "Iter %d: \tk_eff = %1.6f \tCMFD k_eff = %1.6f"
		     "\tres = %1.3E", _num_iterations, _k_eff, cmfd_k_eff, residual);
	else
	  log_printf(NORMAL, "Iter %d: \trho = $%1.4f"
		     "\tres = %1.3E", _num_iterations, getReactivity(), residual);

	/* Check for convergence of the fission source distribution */
	if (i >= 1 && residual < _source_convergence_thresh) {	    
	  _cmfd->getMesh()->computeXS();
	  _cmfd->getMesh()->computeDs(1.0);
	  _timer->stopTimer();
	  _timer->recordSplit("MOC solve");
	  return _k_eff;
	}
    }

    log_printf(WARNING, "Unable to converge the source after %d iterations",
	       max_iterations);

    return _k_eff;
}


/**
 * @brief Computes the total source (fission and scattering) in each flat 
 *        source region.
 * @details This method computes the total source in each region based on
 *          this iteration's current approximation to the scalar flux. A
 *          residual for the source with respect to the source compute on
 *          the previous iteration is computed and returned. The residual
 *          is determined as follows:
 *          /f$ res = \sqrt{\frac{\displaystyle\sum \displaystyle\sum 
 *                    \left(\frac{Q^i - Q^{i-1}{Q^i}\right)^2}{\# FSRs}}} \f$
 *
 * @return the residual between this source and the previous source
 */
FP_PRECISION ThreadPrivateSolverTransient::computeFSRSources() {
    
    int tid;
    FP_PRECISION scatter_source;
    FP_PRECISION fission_source;
    FP_PRECISION delayed_source;
    FP_PRECISION method_source;
    double* nu_sigma_f;
    double* sigma_s;
    double* sigma_t;
    double* chi;
    Material* material;
    Mesh* mesh = _geometry->getMesh();
    Mesh* geom_mesh = _geometry->getGeomMesh();
    int r;
    double fis = 0.0;
    double delay = 0.0;
    double method = 0.0;
    double scat = 0.0;
    double total_source = 0.0;

    FP_PRECISION source_residual = 0.0;
    
    /* For all regions, find the source */
    std::vector<int>::iterator iter;
    double* fine_shape = geom_mesh->getFluxes(PREVIOUS_CONV);
    double* coarse_flux_new = mesh->getFluxes(PREVIOUS);
    double* coarse_flux_old = mesh->getFluxes(PREVIOUS_CONV);
    double* volumes = mesh->getVolumes();
    double* velocity = mesh->getVelocity();
    double* frequency_new = mesh->getFrequencies(CURRENT);
    double* frequency_old = mesh->getFrequencies(PREVIOUS);
    double old_flux;
    double dt = mesh->getDtMOC();
    TimeStepper* ts = mesh->getTimeStepper();
    
    /* loop over coarse mesh cells */
    for (int i = 0; i < mesh->getNumCells(); i++){
	
	/* loop over fsrs within cell */
	for (iter = mesh->getCellFSRs()->at(i).begin(); iter != mesh->getCellFSRs()->at(i).end(); ++iter){
	    
	    r = *iter;

	    material = _FSR_materials[r];
	    nu_sigma_f = material->getNuSigmaF();
	    chi = material->getChi();
	    sigma_s = material->getSigmaS();
	    sigma_t = material->getSigmaT();
	    
	    /* Compute fission source for each group */
	    for (int e=0; e < _num_groups; e++)
		_fission_sources(r,e) = _scalar_flux(r,e) * nu_sigma_f[e];
	    
	    fission_source = pairwise_sum<FP_PRECISION>(&_fission_sources(r,0), 
							_num_groups);
	    
	    /* Compute total scattering source for group G */
	    for (int G=0; G < _num_groups; G++) {
		scatter_source = 0;
		_source_residuals(r,G) = 0.0;
		
		for (int g=0; g < _num_groups; g++)
		    _scatter_sources(r,g) = sigma_s[G*_num_groups+g]*_scalar_flux(r,g);
		
		scatter_source = pairwise_sum<FP_PRECISION>(&_scatter_sources(r,0),
							    _num_groups);
	
		/* delayed neutron precursor source */
		if (mesh->getInitialState() == false){
		    delayed_source = 0.0;
		    if (material->isFissionable()){
			for (int dg=0; dg < mesh->getNumDelayGroups(); dg++)
			    delayed_source += mesh->getLambda()[dg] * 
				static_cast<FunctionalMaterial*>(material)->getPrecConc(CURRENT, dg);
		    }

		    /* method source */
		    method_source = 0.0;
		    old_flux = fine_shape[r*_num_groups+G] * coarse_flux_old[i*_num_groups+G] 
			* volumes[i] / velocity[G];

		    if (mesh->getTransientType() == THETA){
		      method_source = 1.0/(velocity[G] * dt) * 
			(old_flux - _scalar_flux(r,G));	
		    }
		    else if (mesh->getTransientType() == MAF){
			method_source = 1.0/(velocity[G]) * 
			    (old_flux * coarse_flux_new[i*_num_groups+G] 
			    / coarse_flux_old[i*_num_groups+G] / dt 
			     - _scalar_flux(r,G) * (1.0 / dt + frequency_new[i*_num_groups+G]));
		    }
		    
		    /* Set the total source for region r in group G */
		    _source(r,G) = ((1.0 - mesh->getBetaSum()) * (1.0 / mesh->getKeff0()) 
				    * fission_source * chi[G] + scatter_source + chi[G] 
				    * delayed_source + method_source) * ONE_OVER_FOUR_PI;	

		    scat += scatter_source * ONE_OVER_FOUR_PI;
		    delay += chi[G] * delayed_source * ONE_OVER_FOUR_PI;
		    method += chi[G] * method_source * ONE_OVER_FOUR_PI;
		    
		    if (mesh->getTransientType() == ADIABATIC){
		        fis += (1.0 - mesh->getBetaSum()) * (1.0 / _k_eff) 
			  * fission_source * chi[G] * ONE_OVER_FOUR_PI;
		    }
		    else{
		        fis += (1.0 - mesh->getBetaSum()) * (1.0 / mesh->getKeff0()) 
			  * fission_source * chi[G] * ONE_OVER_FOUR_PI;
		    }
		}
		else{
		    /* Set the total source for region r in group G */
		    _source(r,G) = ((1.0 / _k_eff) * fission_source * chi[G]
				    + scatter_source) * ONE_OVER_FOUR_PI;		 

		    total_source += _source(r,G);
		}
		
		_reduced_source(r,G) = _source(r,G) / sigma_t[G];
		
		/* Compute the norm of residual of the source in the region, group */
		if (fabs(_source(r,G)) > 1E-10)
		    _source_residuals(r,G) = pow((_source(r,G) - _old_source(r,G)) 
						 / _source(r,G), 2);
		
		/* Update the old source */
		_old_source(r,G) = _source(r,G);		
	    }
	}
    }
	
    /* Sum up the residuals from each group and in each region */
    source_residual = pairwise_sum<FP_PRECISION>(_source_residuals, 
						 _num_FSRs*_num_groups);
    source_residual = sqrt(source_residual / (_num_groups*_num_FSRs));
    
    return source_residual;
}


/**
 * @brief Allocates memory for track boundary angular fluxes and leakages
 *        flat source region scalar fluxes.
 * @details Deletes memory for old flux arrays if they were allocated from
 *          previous simulation.
 */
void ThreadPrivateSolverTransient::initializeFluxArrays() {

    CPUSolver::initializeFluxArrays();
   
    /* Delete old flux arrays if they exist */
    if (_thread_flux != NULL)
        delete [] _thread_flux;

    /* Delete old current arrays if they exist */
    if (_thread_currents != NULL)
        delete [] _thread_currents;

    int size;

    /* Allocate memory for the flux and leakage arrays */
    try{
	/* Allocate a thread local array of FSR scalar fluxes */
	size = _num_threads * _num_FSRs * _num_groups * sizeof(FP_PRECISION);
	_thread_flux = new FP_PRECISION[size];

	/* Allocate a thread local array of mesh cell surface currents */
	if (_cmfd->getMesh()->getCmfdOn()){ 
	  size = _num_threads * _num_mesh_cells * 8 * _num_groups * sizeof(double);
	  _thread_currents = new double[size];
	}
    }
    catch(std::exception &e) {
        log_printf(ERROR, "Could not allocate memory for the solver's fluxes. "
		   "Backtrace:%s", e.what());
    }
}


/**
 * @brief Set the scalar flux for each energy group inside each flat source 
 *        region to a constant value.
 * @details This method also flattens the thread private flat source region
 *          scalar flux array.
 * @param value the value to assign to each flat source region flux
 */
void ThreadPrivateSolverTransient::flattenFSRFluxes(FP_PRECISION value) {

    CPUSolver::flattenFSRFluxes(value);

    /* Flatten the thread private flat source region scalar flux array */
    #pragma omp parallel for schedule(guided)
    for (int tid=0; tid < _num_threads; tid++) {
        for (int r=0; r < _num_FSRs; r++) {
	    for (int e=0; e < _num_groups; e++) {
	        _thread_flux(tid,r,e) = 0.0;
	    }
        }
    }

    return;
}


/**
 * @brief Set the surface currents for each energy group inside each
 * 			mesh cell to zero
 */
void ThreadPrivateSolverTransient::zeroSurfaceCurrents() {

	CPUSolver::zeroSurfaceCurrents();

    #pragma omp parallel for schedule(guided)
	for (int tid=0; tid < _num_threads; tid++){
		for (int r=0; r < _num_mesh_cells; r++) {
			for (int s=0; s < 8; s++) {
				for (int e=0; e < _num_groups; e++)
					_thread_currents(tid,r*8+s,e) = 0.0;
			}
		}
	}

    return;
}


/**
 * @brief This method performs one transport sweep of all azimuthal angles, 
 *        tracks, segments, polar angles and energy groups.
 * @details The method integrates the flux along each track and updates the 
 *          boundary fluxes for the corresponding output track, while updating 
 *          the scalar flux in each flat source region
 */
void ThreadPrivateSolverTransient::transportSweep() {

    int tid;
    int fsr_id;
    Track* curr_track;
    int azim_index;
    int num_segments;
    segment* curr_segment;    
    segment* segments;
    FP_PRECISION* track_flux;

    log_printf(DEBUG, "Transport sweep with %d OpenMP threads", _num_threads);

    /* Initialize flux in each region to zero */
    flattenFSRFluxes(0.0);

    if (_cmfd->getMesh()->getCmfdOn()) 
      zeroSurfaceCurrents();

    /* Loop over azimuthal angle halfspaces */
    for (int i=0; i < 2; i++) {

        /* Compute the minimum and maximum track IDs corresponding to 
         * this azimuthal angular halfspace */
        int min = i * (_tot_num_tracks / 2);
	int max = (i + 1) * (_tot_num_tracks / 2);
	
	/* Loop over each thread within this azimuthal angle halfspace */
	#pragma omp parallel for private(tid, fsr_id, curr_track, \
	  azim_index, num_segments, segments, curr_segment, \
	  track_flux) schedule(guided)
	for (int track_id=min; track_id < max; track_id++) {

	    tid = omp_get_thread_num();

	    /* Initialize local pointers to important data structures */	
	    curr_track = _tracks[track_id];
	    azim_index = curr_track->getAzimAngleIndex();
	    num_segments = curr_track->getNumSegments();
	    segments = curr_track->getSegments();
	    track_flux = &_boundary_flux(track_id,0,0,0);

	    /* Loop over each segment in forward direction */
	    for (int s=0; s < num_segments; s++) {
	        curr_segment = &segments[s];
		fsr_id = curr_segment->_region_id;
		scalarFluxTally(curr_segment, azim_index,track_flux, 
				&_thread_flux(tid,fsr_id,0),true);
	    }

	    /* Transfer flux to outgoing track */
	    transferBoundaryFlux(track_id, azim_index, true, track_flux);
	    
	    /* Loop over each segment in reverse direction */
	    track_flux += _polar_times_groups;
	    
	    for (int s=num_segments-1; s > -1; s--) {
	        curr_segment = &segments[s];
		fsr_id = curr_segment->_region_id;
		scalarFluxTally(curr_segment, azim_index, track_flux, 
	                        &_thread_flux(tid,fsr_id,0),false);
	    }
	    
	    /* Transfer flux to outgoing track */
	    transferBoundaryFlux(track_id, azim_index, false, track_flux);
	}
    }

    reduceThreadScalarFluxes();

    if (_cmfd->getMesh()->getCmfdOn())
        reduceThreadSurfaceCurrents();
    
    return;
}


/**
 * @brief Computes the contribution to the flat source region scalar flux
 *        from a single track segment.
 * @details This method integrates the angular flux for a track segment across
 *        energy groups and polar angles, and tallies it into the flat
 *        source region scalar flux, and updates the track's angular flux.
 * @param curr_segment a pointer to the segment of interest
 * @param track_flux a pointer to the track's angular flux
 * @param fsr_flux a pointer to the temporary flat source region flux buffer
 */
void ThreadPrivateSolverTransient::scalarFluxTally(segment* curr_segment, 
	                                           int azim_index,
						   FP_PRECISION* track_flux,
						   FP_PRECISION* fsr_flux,
						   bool fwd){

    int tid = omp_get_thread_num();
    int fsr_id = curr_segment->_region_id;
    FP_PRECISION length = curr_segment->_length;
    double* sigma_t = curr_segment->_material->getSigmaT();

    /* The average flux along this segment in the flat source region */
    FP_PRECISION deltapsi;
    FP_PRECISION exponential;

    /* Loop over energy groups */
    for (int e=0; e < _num_groups; e++) {

	/* Loop over polar angles */
        for (int p=0; p < _num_polar; p++){
            exponential = computeExponential(sigma_t[e], length, p);
            deltapsi = (track_flux(p,e) - _reduced_source(fsr_id,e)) * exponential;
	    fsr_flux[e] += deltapsi * _polar_weights(azim_index, p);
	    track_flux(p,e) -= deltapsi;
	}
    }

    if (_cmfd->getMesh()->getCmfdOn()){
    	if (curr_segment->_mesh_surface_fwd != -1 && fwd){

    		/* set polar angle * energy group to 0 */
    		int pe = 0;

    		/* loop over energy groups */
    		for (int e = 0; e < _num_groups; e++) {

    			/* loop over polar angles */
    			for (int p = 0; p < _num_polar; p++){

    				/* increment current (polar and azimuthal weighted flux, group)*/
			    _thread_currents(tid,curr_segment->_mesh_surface_fwd,e) += track_flux(p,e)*_polar_weights(azim_index, p)/2.0;

    				pe++;
    			}
    		}
    	}
    	else if (curr_segment->_mesh_surface_bwd != -1 && !fwd){

    		/* set polar angle * energy group to 0 */
    		int pe = 0;

    		/* loop over energy groups */
    		for (int e = 0; e < _num_groups; e++) {

    			/* loop over polar angles */
    			for (int p = 0; p < _num_polar; p++){

    				/* increment current (polar and azimuthal weighted flux, group)*/
			    _thread_currents(tid,curr_segment->_mesh_surface_bwd,e) += track_flux(p,e)*_polar_weights(azim_index, p)/2.0;

    				pe++;
    			}
    		}
    	}
    }

    return;
}


/**
 * @brief Reduces the flat source region scalar fluxes from private thread 
 *        array to a global array.
 */
void ThreadPrivateSolverTransient::reduceThreadScalarFluxes() {

    for (int tid=0; tid < _num_threads; tid++) {
        for (int r=0; r < _num_FSRs; r++) {
            for (int e=0; e < _num_groups; e++) {
                _scalar_flux(r,e) += _thread_flux(tid,r,e);
            }
        }
    }

    return;
}


/**
 * @brief Reduces the surface currents from private thread
 *        array to a global array.
 */
void ThreadPrivateSolverTransient::reduceThreadSurfaceCurrents() {

	for (int tid=0; tid < _num_threads; tid++){
		for (int r=0; r < _num_mesh_cells; r++) {
			for (int s=0; s < 8; s++) {
				for (int e=0; e < _num_groups; e++)
					_surface_currents(r*8+s,e) += _thread_currents(tid,r*8+s,e);
			}
		}
	}

    return;
}


void ThreadPrivateSolverTransient::scaleTrackFlux(double scale_val){

    #pragma omp parallel for schedule(guided)
    for (int t=0; t < _tot_num_tracks; t++) {
        for (int d=0; d < 2; d++) {
            for (int p=0; p < _num_polar; p++) {
	        for (int e=0; e < _num_groups; e++) {
		    _boundary_flux(t,d,p,e) *= scale_val;
	        }  
	    }
        }
    }
}


double ThreadPrivateSolverTransient::getReactivity(){

  Mesh* mesh = _geometry->getMesh();
  double reactivity = (_k_eff - mesh->getKeff0()) / mesh->getKeff0() / mesh->getBetaSum();

  return reactivity;
}


int ThreadPrivateSolverTransient::getNumIters(){
  return _num_iterations;
}
