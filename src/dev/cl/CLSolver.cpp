#include "CLSolver.h"

/**
 * @brief Constructor initializes array pointers for tracks and materials.
 * @details The constructor retrieves the number of energy groups and flat
 *          source regions and azimuthal angles from the geometry and track
 *          generator.
 * @param geometry an optional pointer to the geometry
 * @param track_generator an optional pointer to the track generator
 */
 
CLSolver::CLSolver(Geometry* geometry, TrackGenerator* track_generator) :
    Solver(geometry, track_generator) {

    /**************************************************************************/
    /*                        Host data initialization                        */
    /**************************************************************************/

    /* The default number of threadblocks and threads per threadblock */
     _B = 64;
    _T = 64;


    /**************************************************************************/
    /*                       Device data initialization                       */
    /**************************************************************************/

    _materials = NULL;
    _dev_tracks = NULL;

    _tot_absorption = NULL;
    _tot_fission = NULL;
    _leakage = NULL;

    if (track_generator != NULL)
        setTrackGenerator(track_generator);

    if (geometry != NULL)
        setGeometry(geometry);

    if(!machineContainsCL()) {
        log_printf(ERROR, "Unable to attach opencl solver since no opencl
                device is attached to the machine");
    } else {
        inst = clInst;
    }
}


/**
 * @brief Solver destructor frees all memory on the device, including arrays
 *        for the fluxes and sources.
 */
CLSolver::~CLSolver() {

    if (_FSR_volumes != NULL) {
        cudaFree(_FSR_volumes);
	_FSR_volumes = NULL;
    }

    if (_FSR_materials != NULL) {
        cudaFree(_FSR_materials);
	_FSR_materials = NULL;
    }

    if (_materials != NULL) {
        cudaFree(_materials);
	_materials = NULL;
    }

    if (_dev_tracks != NULL) {
        cudaFree(_dev_tracks);
	_dev_tracks = NULL;
    }

    if (_boundary_flux != NULL) {
        cudaFree(_boundary_flux);
	_boundary_flux = NULL;
    }

    if (_scalar_flux != NULL) {
        cudaFree(_scalar_flux);
	_scalar_flux = NULL;
    }

    if (_source != NULL) {
        cudaFree(_source);
	_source = NULL;
    }

    if (_old_source != NULL) {
        cudaFree(_old_source);
	_old_source = NULL;
    }

    if (_reduced_source != NULL) {
        cudaFree(_reduced_source);
	_reduced_source = NULL;
    }

    if (_fission_sources != NULL) {
        _fission_sources_vec.clear();
	_fission_sources = NULL;
    }

    if (_tot_absorption != NULL) {
        _tot_absorption_vec.clear();
	_tot_absorption = NULL;
    }

    if (_tot_fission != NULL) {
        _tot_fission_vec.clear();
	_tot_fission = NULL;
    }

    if (_source_residuals != NULL) {
        _source_residuals_vec.clear();
	_source_residuals = NULL;
    }

    if (_leakage != NULL) {
        _leakage_vec.clear();
	_leakage = NULL;
    }

    if (_prefactor_array != NULL) {
        cudaFree(_prefactor_array);
	_prefactor_array = NULL;
    }
}


/**
 * @brief Returns the scalar flux for some energy group for a flat source region
 * @param fsr_id the ID for the FSR of interest
 * @param energy_group the energy group of interest
 */
FP_PRECISION CLSolver::getFSRScalarFlux(int fsr_id, int energy_group) {

    /* Error checking */
    if (fsr_id >= _num_FSRs)
        log_printf(ERROR, "Unable to return a scalar flux for FSR id = %d "
		 "in enery group %d since the solver only contains FSR with "
		   "IDs greater than or equal to %d", 
		   fsr_id, energy_group, _num_FSRs-1);

    if (fsr_id < 0)
        log_printf(ERROR, "Unable to return a scalar flux for FSR id = %d "
		  "in energy group %d since FSRs do not have negative IDs", 
		  fsr_id, energy_group);

    if (energy_group-1 >= _num_groups)
        log_printf(ERROR, "Unable to return a scalar flux for FSR id = %d "
		   "in energy group %d since the solver only has %d energy "
		   "groups", fsr_id, energy_group, _num_groups);

    if (energy_group <= 0)
        log_printf(ERROR, "Unable to return a scalar flux for FSR id = %d "
		 "in energy group %d since energy groups are greater than 1",
		 fsr_id, energy_group);

    /* Copy the scalar flux for this FSR and energy group from the device */
    FP_PRECISION fsr_scalar_flux;
    int flux_index = fsr_id * _num_groups + energy_group - 1;
    cudaMemcpy((void*)&fsr_scalar_flux, (void*)&_scalar_flux[flux_index], 
	       sizeof(FP_PRECISION), cudaMemcpyDeviceToHost);

    return fsr_scalar_flux;
}


/**
 * @brief Return an array indexed by flat source region IDs and energy groups 
 *        which contains the corresponding fluxes for each flat source region.
 * @return an array of flat source region scalar fluxes
 */
FP_PRECISION* CLSolver::getFSRScalarFluxes() {

    if (_scalar_flux == NULL)
        log_printf(ERROR, "Unable to returns the device solver's scalar flux "
		   "array since it has not yet been allocated in memory");

    printf("Retrieving FSR scalar fluxes from GPU\n");

    /* Copy the scalar flux for all FSRs from the device */
    FP_PRECISION* fsr_scalar_fluxes = new FP_PRECISION[_num_FSRs * _num_groups];
    cudaMemcpy((void*)fsr_scalar_fluxes, (void*)_scalar_flux,
	       _num_FSRs * _num_groups * sizeof(FP_PRECISION),
	       cudaMemcpyDeviceToHost);

    printf("Finished retrieving, returning\n");

    return fsr_scalar_fluxes;
}


/**
 * @brief Sets the number of threadblocks (>0) for device kernels
 * @param num_blocks the number of threadblocks
 */
void CLSolver::setNumThreadBlocks(int num_blocks) {

    if (num_blocks < 0)
        log_printf(ERROR, "Unable to set the number of threadblocks to %d since "
		   "it is a negative number", num_blocks);

    _B = num_blocks;
}


/**
 * @brief Sets the number of threads per block (>0) for device kernels
 * @param num_threads the number of threads per block
 */
void CLSolver::setNumThreadsPerBlock(int num_threads) {

    if (num_threads < 0)
        log_printf(ERROR, "Unable to set the number of threads per block to %d "
		   "since it is a negative number", num_threads);

    _T = num_threads;
}


/**
 * @brief Sets the geometry for the solver.
 * @details The geometry must already have initialized flat source region maps
 *          and segmentized the trackgenerator's tracks.
 * @param geometry a pointer to a geometry
 */
void CLSolver::setGeometry(Geometry* geometry) {
    Solver::setGeometry(geometry);
    initializeMaterials();

    /* Copy the number of energy groups to constant memory on the GPU */
    cudaMemcpyToSymbol(num_groups, (void*)&_num_groups, sizeof(int), 0,
		       cudaMemcpyHostToDevice);
}


/**
 * @brief Sets the track generator with characteristic tracks for the solver.
 * @details The track generator must already have generated tracks and have
 *          segmentized them across the geometry.
 * @param track_generator a pointer to a trackgenerator
 */
void CLSolver::setTrackGenerator(TrackGenerator* track_generator) {
    Solver::setTrackGenerator(track_generator);
    initializeTracks();
}


/**
 * @brief Creates a polar quadrature object for the solver.
 */
void CLSolver::initializePolarQuadrature() {

    log_printf(INFO, "Initializing polar quadrature on the GPU...");

    /* Deletes the old quadrature if one existed */
    if (_quad != NULL)
        delete _quad;

    _quad = new Quadrature(_quadrature_type, _num_polar);
    _polar_times_groups = _num_groups * _num_polar;

    /* Copy the number of polar angles to constant memory on the GPU */
    cudaMemcpyToSymbol(num_polar, (void*)&_num_polar, sizeof(int), 0,
		       cudaMemcpyHostToDevice);

    /* Copy twice the number of polar angles to constant memory on the GPU */
    cudaMemcpyToSymbol(two_times_num_polar, (void*)&_two_times_num_polar, 
		       sizeof(int), 0, cudaMemcpyHostToDevice);

    /* Copy the number of polar angles times energy groups to constant memory 
     * on the GPU */
    cudaMemcpyToSymbol(polar_times_groups, (void*)&_polar_times_groups, 
		       sizeof(int), 0, cudaMemcpyHostToDevice);

    /* Compute polar times azimuthal angle weights */
    if (_polar_weights != NULL)
        delete [] _polar_weights;

    _polar_weights =
        (FP_PRECISION*)malloc(_num_polar * _num_azim * sizeof(FP_PRECISION));

    FP_PRECISION* multiples = _quad->getMultiples();
    double* azim_weights = _track_generator->getAzimWeights();

    for (int i=0; i < _num_azim; i++) {
        for (int j=0; j < _num_polar; j++)
	    _polar_weights[i*_num_polar+j] = azim_weights[i]*multiples[j]*FOUR_PI;
    }

    /* Copy the polar weights to constant memory on the GPU */
    cudaMemcpyToSymbol(polar_weights, (void*)_polar_weights,
		       _num_polar * _num_azim * sizeof(FP_PRECISION),
		       0, cudaMemcpyHostToDevice);
}


/**
 * @brief Initializes each of the flat source region objects inside the solver's
 *        array of flatsourceregions. 
 * @details This method assigns each flat source region a unique, monotonically
 *          increasing ID, sets the material for each flat source region, and 
 *          assigns a volume based on the cumulative length of all of the 
 *          segments inside the flat source region.
 */
void CLSolver::initializeFSRs() {

    log_printf(INFO, "Initializing FSRs on the GPU...");

    /* Delete old FSRs array if it exists */
    if (_FSR_volumes != NULL)
        cudaFree(_FSR_volumes);

    if (_FSR_materials != NULL)
        cudaFree(_FSR_materials);

    /* Allocate memory for all tracks and track offset indices on the device */
    try{

        /* Allocate memory on device for FSR volumes and material uids */
        cudaMalloc((void**)&_FSR_volumes, _num_FSRs * sizeof(FP_PRECISION));
        cudaMalloc((void**)&_FSR_materials, _num_FSRs * sizeof(int));
	
	/* Create a temporary FSR array to populate and then copy to device */
	FP_PRECISION* temp_FSR_volumes = new FP_PRECISION[_num_FSRs];

	/* Get the array indexed by FSR IDs with material ID values */
	int* FSRs_to_materials = _geometry->getFSRtoMaterialMap();

	/* Initialize each FSRs volume to 0 to avoid NaNs */
	memset(temp_FSR_volumes, FP_PRECISION(0.), 
	       _num_FSRs*sizeof(FP_PRECISION));

	Track* track;
	int num_segments;
	segment* curr_segment;
	segment* segments;
	FP_PRECISION volume;

	double* azim_weights = _track_generator->getAzimWeights();

	/* Set each FSR's volume by accumulating the total length of all
	   tracks inside the FSR. Iterate over azimuthal angle, track, segment*/
	for (int i=0; i < _num_azim; i++) {
	    for (int j=0; j < _num_tracks[i]; j++) {

	        track = &_track_generator->getTracks()[i][j];
		num_segments = track->getNumSegments();
		segments = track->getSegments();

		/* Iterate over the track's segments to update FSR volumes */
		for (int s = 0; s < num_segments; s++) {
		    curr_segment = &segments[s];
		    volume = curr_segment->_length * azim_weights[i];
		    temp_FSR_volumes[curr_segment->_region_id] += volume;
		}
	    }
	}

	/* Copy the temporary array of FSRs to the device */
	cudaMemcpy((void*)_FSR_volumes, (void*)temp_FSR_volumes,
		   _num_FSRs * sizeof(FP_PRECISION), cudaMemcpyHostToDevice);
	cudaMemcpy((void*)_FSR_materials, (void*)FSRs_to_materials,
		   _num_FSRs * sizeof(int), cudaMemcpyHostToDevice);

	/* Copy the number of flat source regions into constant memory on 
	 * the GPU */
	cudaMemcpyToSymbol(num_FSRs, (void*)&_num_FSRs, sizeof(int), 0,
		       cudaMemcpyHostToDevice);

	/* Free the temporary array of FSRs on the host */
	free(temp_FSR_volumes);
    }
    catch(std::exception &e) {
        log_printf(ERROR, "Could not allocate memory for the solver's flat "
		   "source regions on the device. Backtrace:%s", e.what());
    }

    initializeThrustVectors();
}


/**
 * @brief Allocates data on the GPU for all materials data.
 */
void CLSolver::initializeMaterials() {

    log_printf(INFO, "Initializing materials on the GPU...");

    /* Delete old materials array if it exists */
    if (_materials != NULL)
        cudaFree(_materials);

    /* Allocate memory for all tracks and track offset indices on the device */
    try{

	std::map<int, Material*> host_materials=_geometry->getMaterials();
	std::map<int, Material*>::iterator iter;

        /* Iterate through all materials and clone them on the device */
        cudaMalloc((void**)&_materials, _num_materials * sizeof(dev_material));
	for (iter=host_materials.begin(); iter != host_materials.end(); ++iter)
	    cloneMaterialOnGPU(iter->second, &_materials[iter->second->getUid()]);
    }
    catch(std::exception &e) {
        log_printf(ERROR, "Could not allocate memory for the device solver's "
		   "materials. Backtrace:%s", e.what());
    }
}


/**
 * @brief Allocates memory on the GPU for all tracks in the simulation.
 */
void CLSolver::initializeTracks() {

    log_printf(INFO, "Initializing tracks on the GPU...");

    /* Delete old tracks array if it exists */
    if (_dev_tracks != NULL)
        cudaFree(_dev_tracks);

    /* Allocate memory for all tracks and track offset indices on the device */
    try{

        /* Allocate array of tracks */
    	cudaMalloc((void**)&_dev_tracks, _tot_num_tracks * sizeof(dev_track));

        /* Iterate through all tracks and clone them on the device */
	int index;

	for (int i=0; i < _tot_num_tracks; i++) {

	    cloneTrackOnGPU(_tracks[i], &_dev_tracks[i]);

	    /* Make track reflective */
	    index = computeScalarTrackIndex(_tracks[i]->getTrackInI(),
		        		       _tracks[i]->getTrackInJ());
	    cudaMemcpy((void*)&_dev_tracks[i]._track_in,
		   (void*)&index, sizeof(int), cudaMemcpyHostToDevice);

	    index = computeScalarTrackIndex(_tracks[i]->getTrackOutI(),
		        		       _tracks[i]->getTrackOutJ());
	    cudaMemcpy((void*)&_dev_tracks[i]._track_out,
		   (void*)&index, sizeof(int), cudaMemcpyHostToDevice);
	}

	/* Copy the array of number of tracks for each azimuthal angles into 
	 * constant memory on GPU */
	cudaMemcpyToSymbol(num_tracks, (void*)_num_tracks, 
			   _num_azim * sizeof(int), 0, cudaMemcpyHostToDevice);
    
	/* Copy the total number of tracks into constant memory on GPU */
	cudaMemcpyToSymbol(tot_num_tracks, (void*)&_tot_num_tracks,
			   sizeof(int), 0, cudaMemcpyHostToDevice);

	/* Copy the number of azimuthal angles into constant memory on GPU */
	cudaMemcpyToSymbol(num_azim, (void*)&_num_azim, sizeof(int), 0, 
			   cudaMemcpyHostToDevice);
	
	/* Copy the array of number of tracks for each azimuthal angles into 
	 * constant memory on GPU */
	cudaMemcpyToSymbol(num_tracks, (void*)_num_tracks, 
			   _num_azim * sizeof(int), 0, cudaMemcpyHostToDevice);
	
	/* Copy the total number of tracks into constant memory on GPU */
	cudaMemcpyToSymbol(tot_num_tracks, (void*)&_tot_num_tracks,
			   sizeof(int), 0, cudaMemcpyHostToDevice);
    }

    catch(std::exception &e) {
        log_printf(ERROR, "Could not allocate memory for the solver's tracks "
		   "on the device. Backtrace:%s", e.what());
    }
}


/**
 * @brief Allocates memory for track boundary angular fluxes and leakages
 *        flat source region scalar fluxes.
 * @details Deletes memory for old flux arrays if they were allocated from
 *          previous simulation.
 */
void CLSolver::initializeFluxArrays() {

    log_printf(INFO, "Initializing flux arrays on the GPU...");

    /* Delete old flux arrays if they exist */
    if (_boundary_flux != NULL)
        cudaFree(_boundary_flux);

    if (_scalar_flux != NULL)
        cudaFree(_scalar_flux);

    /* Allocate memory for all flux arrays on the device */
    try{
        cudaMalloc((void**)&_boundary_flux,
		   2*_tot_num_tracks * _polar_times_groups*sizeof(FP_PRECISION));
        cudaMalloc((void**)&_scalar_flux, 
		   _num_FSRs * _num_groups * sizeof(FP_PRECISION));
    }
    catch(std::exception &e) {
        log_printf(ERROR, "Could not allocate memory for the solver's fluxes "
		   "on the device. Backtrace:%s", e.what());
    }
}


/**
 * @brief Allocates memory for flat source region source arrays.
 * @details Deletes memory for old source arrays if they were allocated from
 *          previous simulation.
 */
void CLSolver::initializeSourceArrays() {

    log_printf(INFO, "Initializing source arrays on the GPU...");

    /* Delete old sources arrays if they exist */
    if (_source != NULL)
        cudaFree(_source);

    if (_old_source != NULL)
        cudaFree(_old_source);

    if (_reduced_source != NULL)
        cudaFree(_reduced_source);

    /* Allocate memory for all source arrays on the device */
    try{

        cudaMalloc((void**)&_source, 
		   _num_FSRs * _num_groups * sizeof(FP_PRECISION));

	cudaMalloc((void**)&_old_source,
		   _num_FSRs * _num_groups * sizeof(FP_PRECISION));

	cudaMalloc((void**)&_reduced_source,
		   _num_FSRs * _num_groups * sizeof(FP_PRECISION));
    }
    catch(std::exception &e) {
        log_printf(ERROR, "Could not allocate memory for the solver's flat "
		   "source region sources array on the device. "
		   "Backtrace:%s", e.what());
    }
}


/**
 * @brief Initialize Thrust vectors for the fission and absorption rates,
 *        source residuals, leakage and fission sources.
 */
void CLSolver::initializeThrustVectors() {

    log_printf(INFO, "Initializing thrust vectors on the GPU...");

    /* Delete old vectors if they exist */
    if (_fission_sources != NULL) {
        _fission_sources = NULL;
        _fission_sources_vec.clear();
    }

    if (_tot_absorption != NULL) {
        _tot_absorption = NULL;
        _tot_absorption_vec.clear();
    }

    if (_tot_fission != NULL) {
        _tot_fission = NULL;
        _tot_fission_vec.clear();
    }

    if (_source_residuals != NULL) {
        _source_residuals = NULL;
        _source_residuals_vec.clear();
    }

    if (_leakage != NULL) {
        _leakage = NULL;
        _leakage_vec.clear();
    }


    /* Allocate memory for fission, absorption and source vectors on device */
    try{
        /* Allocate fission source array on device */
        _fission_sources_vec.resize(_B * _T);
	_fission_sources = thrust::raw_pointer_cast(&_fission_sources_vec[0]);
      
	/* Allocate total absorption reaction rate array on device */
	_tot_absorption_vec.resize(_B * _T);
	_tot_absorption = thrust::raw_pointer_cast(&_tot_absorption_vec[0]);

	/* Allocate fission reaction rate array on device */
	_tot_fission_vec.resize(_B * _T);
	_tot_fission = thrust::raw_pointer_cast(&_tot_fission_vec[0]);

	/* Allocate source residual array on device */
	_source_residuals_vec.resize(_B * _T);
	_source_residuals = thrust::raw_pointer_cast(&_source_residuals_vec[0]);

	/* Allocate leakage array on device */
	_leakage_vec.resize(_B * _T);
	_leakage = thrust::raw_pointer_cast(&_leakage_vec[0]);
    }
    catch(std::exception &e) {
        log_printf(ERROR, "Could not allocate memory for the solver's "
		   "Thrust vectors.Backtrace:%s", e.what());
    }
}


/**
 * @brief This method computes the index for the jth track at azimuthal angle i.
 * @details This method is necessary since the array of tracks on the device 
 *          is a 1D array which needs a one-to-one mapping from the 2D jagged 
 *          array of tracks on the host.
 * @param i azimuthal angle number
 * @param j the jth track at angle i
 * @return an index into the device track array
 */
int CLSolver::computeScalarTrackIndex(int i, int j) {

    int index =0;
    int p = 0;

    /* Iterate over each azimuthal angle and increment index by the number of
       tracks at each angle */
    while (p < i) {
        index += _num_tracks[p];
	p++;
    }

    /* Update index for this track since it is the jth track at angle i */
    index += j;
    
    return index;
}


/**
 * @brief Pre-computes exponential pre-factors for each segment of each track 
 *        for each polar angle. 
 * @details This method will generate a hashmap which contains values of the 
 *          prefactor for specific segment lengths (the keys into the hashmap).
 */
void CLSolver::precomputePrefactors(){

    log_printf(INFO, "Building exponential prefactor hashtable on device...");

    /* Copy a boolean indicating whether or not to use the linear interpolation
     * table or the exp intrinsic function */
    cudaMemcpyToSymbol(interpolate_exponential, 
		       (void*)&_interpolate_exponential, 
		       sizeof(bool), 0, cudaMemcpyHostToDevice);

    /* Copy the sines of the polar angles which is needed if the user
     * requested the use of the exp intrinsic to evaluate exponentials */
    cudaMemcpyToSymbol(sinthetas, (void*)_quad->getSinThetas(),
		       _num_polar * sizeof(FP_PRECISION), 0, 
		       cudaMemcpyHostToDevice);

    /* Set size of prefactor array */
    int num_array_values = 10 * sqrt(1. / (8. * _source_convergence_thresh));
    _prefactor_spacing = 10. / num_array_values;
    _inverse_prefactor_spacing = 1.0 / _prefactor_spacing;
    _prefactor_array_size = _two_times_num_polar * num_array_values;
    _prefactor_max_index = _prefactor_array_size - _two_times_num_polar - 1;
    
    /* allocate arrays */
    FP_PRECISION* prefactor_array = new FP_PRECISION[_prefactor_array_size];
    
    FP_PRECISION expon;
    FP_PRECISION intercept;
    FP_PRECISION slope;

    /* Create prefactor array */
    for (int i = 0; i < num_array_values; i ++){
        for (int p = 0; p < _num_polar; p++){
	    expon = exp(- (i * _prefactor_spacing) / _quad->getSinTheta(p));
	    slope = - expon / _quad->getSinTheta(p);
	    intercept = expon * (1 + (i * _prefactor_spacing) /
				 _quad->getSinTheta(p));
	    prefactor_array[_two_times_num_polar * i + 2 * p] = slope;
	    prefactor_array[_two_times_num_polar * i + 2 * p + 1] = intercept;
	}
    }

    /* Allocate memory for the prefactor array on the device */
    cudaMalloc((void**)&_prefactor_array, 
	       _prefactor_array_size * sizeof(FP_PRECISION));

    /* Copy prefactor array to the device */
    cudaMemcpy((void*)_prefactor_array, (void*)prefactor_array, 
	       _prefactor_array_size * sizeof(FP_PRECISION),
	       cudaMemcpyHostToDevice);

    /* Copy prefactor array size and spacing to constant memory on the device */
    cudaMemcpyToSymbol(prefactor_spacing, (void*)&_prefactor_spacing, 
		       sizeof(FP_PRECISION), 0, cudaMemcpyHostToDevice);

    cudaMemcpyToSymbol(inverse_prefactor_spacing, 
		       (void*)&_inverse_prefactor_spacing, 
		       sizeof(FP_PRECISION), 0, cudaMemcpyHostToDevice);

    cudaMemcpyToSymbol(prefactor_max_index, (void*)&_prefactor_max_index,
		       sizeof(int), 0, cudaMemcpyHostToDevice);

    free(prefactor_array);

    return;
}


/**
 * @brief Zero each track's boundary fluxes for each energy group and polar
 *        angle in the "forward" and "reverse" directions.
 */
void CLSolver::zeroTrackFluxes() {
    int size = 2 * _tot_num_tracks * _num_polar * _num_groups;
    size *= sizeof(FP_PRECISION);
    cudaMemset(_boundary_flux, 0.0, size);
    return;
}


/**
 * @brief Set the scalar flux for each energy group inside each flat source 
 *        region to a constant value.
 * @param value the value to assign to each flat source region flux
 */
void CLSolver::flattenFSRFluxes(FP_PRECISION value) {
    int size = _num_FSRs * _num_groups * sizeof(FP_PRECISION);
    cudaMemset(_scalar_flux, value, size);
    return;
}


/**
 * @brief Set the source for each energy group inside each flat source region
 *        to a constant value.
 * @param value the value to assign to each flat source region source
 */
void CLSolver::flattenFSRSources(FP_PRECISION value) {
    int size = _num_FSRs * _num_groups * sizeof(FP_PRECISION);
    cudaMemset(_source, value, size);
    cudaMemset(_old_source, value, size);
    return;
}


/**
 * @brief Normalizes all flat source region scalar fluxes and track boundary
 *        angular fluxes to the total fission source (times \f$ \nu \f$).
 */
void CLSolver::normalizeFluxes() {

    int shared_mem = sizeof(FP_PRECISION) * _T;

    computeFissionSourcesOnDevice<<<_B, _T, shared_mem>>>(_FSR_volumes, 
							  _FSR_materials,
							  _materials, 
							  _scalar_flux, 
							  _fission_sources);

    FP_PRECISION norm_factor = 1.0 / thrust::reduce(_fission_sources_vec.begin(),
						    _fission_sources_vec.end());

    normalizeFluxesOnDevice<<<_B, _T>>>(_scalar_flux, _boundary_flux, 
					norm_factor);
}


FP_PRECISION CLSolver::computeFSRSources() {
    
    computeFSRSourcesOnDevice<<<_B, _T>>>(_FSR_materials, _materials, 
					  _scalar_flux, _source, 
					  _old_source, _reduced_source,
					  1.0 / _k_eff, _source_residuals);

    FP_PRECISION residual = thrust::reduce(_source_residuals_vec.begin(), 
					   _source_residuals_vec.end());
    residual = sqrt(residual / _num_FSRs);

    return residual;
}



void CLSolver::transportSweep() {

    int shared_mem = _T * _two_times_num_polar * sizeof(FP_PRECISION);
    int tid_offset, tid_max;

    log_printf(DEBUG, "Transport sweep on device with %d blocks" 
                      " and %d threads", _B, _T);

    /* Initialize leakage to zero */
    thrust::fill(_leakage_vec.begin(), _leakage_vec.end(), 0.0);
    
    /* Initialize flux in each region to zero */
    flattenFSRFluxes(0.0);
    
    /* Sweep the first halfspace of azimuthal angle space */
    tid_offset = 0;
    tid_max = (_tot_num_tracks / 2);

    transportSweepOnDevice<<<_B, _T, shared_mem>>>(_scalar_flux, 
						   _boundary_flux,
						   _reduced_source, _leakage,
						   _materials, _dev_tracks,
						   _prefactor_array, 
						   tid_offset, tid_max);

    /* Sweep the second halfspace of azimuthal angle space */
    tid_offset = tid_max * _num_groups;
    tid_max = _tot_num_tracks;

    transportSweepOnDevice<<<_B, _T, shared_mem>>>(_scalar_flux,
						   _boundary_flux,
						   _reduced_source, _leakage,
						   _materials, _dev_tracks,
						   _prefactor_array,
						   tid_offset, tid_max);
}


/**
 * @brief Add the source term contribution in the transport equation to 
 *        the flat source region scalar flux
 */
void CLSolver::addSourceToScalarFlux() {

    addSourceToScalarFluxOnDevice<<<_B,_T>>>(_scalar_flux, _reduced_source,
					     _FSR_volumes, _FSR_materials,
					     _materials);
}


/**
 * @brief Compute \f$ k_{eff} \f$ from the total fission and absorption rates.
 * @details This method computes the current approximation to the 
 *          multiplication factor on this iteration as follows:
 *          \f$ k_{eff} = \frac{\displaystyle\sum \displaystyle\sum \nu
 *                        \Sigma_f \Phi V}{\displaystyle\sum 
 *                        \displaystyle\sum \Sigma_a \Phi V} \f$
 */
void CLSolver::computeKeff() {

    FP_PRECISION tot_absorption;
    FP_PRECISION tot_fission;
    FP_PRECISION tot_leakage;

    /* Compute the total fission and absorption rates on the device.
     * This kernel stores partial rates in a Thrust vector with as many
     * entries as GPU threads executed by the kernel */
    computeFissionAndAbsorption<<<_B, _T>>>(_FSR_volumes, _FSR_materials,
					    _materials, _scalar_flux,
					    _tot_absorption, _tot_fission);

    cudaDeviceSynchronize();

    /* Compute the total absorption rate by reducing the partial absorption
     * rates compiled in the Thrust vector */
    tot_absorption = thrust::reduce(_tot_absorption_vec.begin(),
				    _tot_absorption_vec.end());

    /* Compute the total fission rate by reducing the partial fission
     * rates compiled in the Thrust vector */
    tot_fission = thrust::reduce(_tot_fission_vec.begin(),
				 _tot_fission_vec.end());

    cudaMemcpy((void*)&tot_fission, (void*)_tot_fission, 
	       _B * _T * sizeof(FP_PRECISION), cudaMemcpyHostToDevice);


    /* Compute the total leakage by reducing the partial leakage
     * rates compiled in the Thrust vector */
    tot_leakage = 0.5 * thrust::reduce(_leakage_vec.begin(),
				       _leakage_vec.end());


    /* Compute the new keff from the fission and absorption rates */
    _k_eff = tot_fission / (tot_absorption + tot_leakage);

    log_printf(DEBUG, "abs = %f, fiss = %f, leak = %f, keff = %f", 
	       tot_absorption, tot_fission, tot_leakage, _k_eff);
}


/**
 * @brief Computes the volume-weighted, energy integrated fission rate in 
 *        each flat source region and stores them in an array indexed by 
 *        flat source region ID.
 * @param fission_rates an array to store the fission rates, passed in as a
 *        numpy array from Python
 * @param num_FSRs the number of FSRs passed in from Python
 */
void CLSolver::computeFSRFissionRates(double* fission_rates, int num_FSRs) {

    log_printf(INFO, "Computing FSR fission rates...");

    /* Allocate memory for the FSR fission rates on the device */
    double* dev_fission_rates;
    cudaMalloc((void**)&dev_fission_rates, _num_FSRs * sizeof(double));

    /* Compute the FSR fission rates on the device */
    computeFSRFissionRatesOnDevice<<<_B,_T>>>(dev_fission_rates, 
					      _FSR_materials,
					      _materials,
					      _scalar_flux);

    /* Copy the fission rate array from the device to the host */
    cudaMemcpy((void*)fission_rates, (void*)dev_fission_rates, 
	       _num_FSRs * sizeof(double), cudaMemcpyDeviceToHost);

    /* Deallocate the memory assigned to store the fission rates on the device */
    cudaFree(dev_fission_rates);

    return;
}
