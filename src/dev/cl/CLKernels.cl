#include "CLSolver.h"


/** The number of azimuthal angles */
__constant int num_azim[1];

/** The number of energy groups */
__constant int num_groups[1];

/** The number of flat source regions */
__constant int num_FSRs[1];

/** The number of polar angles */
__constant int num_polar[1];

/** Twice the number of polar angles */
__constant int two_times_num_polar[1];

/** The number of polar angles times energy groups */
__constant int polar_times_groups[1];

/** An array for the sines of the polar angle in the polar quadrature set */
__constant FP_PRECISION sinthetas[MAX_POLAR_ANGLES];

/** An array of the weights for the polar angles from the quadrature set */
__constant FP_PRECISION polar_weights[MAX_POLAR_ANGLES*MAX_AZIM_ANGLES];

/** A pointer to an array with the number of tracks per azimuthal angle */
__constant int num_tracks[MAX_AZIM_ANGLES/2];

/** The total number of tracks */
__constant int tot_num_tracks[1];

__constant bool interpolate_exponential[1];

/** The maximum index of the exponential prefactor array */
__constant int prefactor_max_index[1];

/** The spacing for the exponential prefactor array */
__constant FP_PRECISION prefactor_spacing[1];

/** The inverse spacing for the exponential prefactor array */
__constant FP_PRECISION inverse_prefactor_spacing[1];



/**
 * @brief Compute the total fission source from all flat source regions.
 * @param FSR_volumes an array of flat source region volumes
 * @param FSR_materials an array of flat source region materials
 * @param materials an array of materials on the device
 * @param scalar_flux the scalar flux in each flat source region
 * @param fission_sources array of fission sources in each flat source region
 */
__kernel void computeFissionSourcesOnDevice(
        __global FP_PRECISION* FSR_volumes,
        __global int* FSR_materials,
        __global dev_material* materials,
        __global FP_PRECISION* scalar_flux,
        __global FP_PRECISION* fission_sources) {

    /* Use a shared memory buffer for each thread's fission source */
    extern __local FP_PRECISION shared_fission_source[];

    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    dev_material* curr_material;
    double* nu_sigma_f;
    FP_PRECISION volume;
    FP_PRECISION source;

    /* Initialize fission source to zero */
    shared_fission_source[threadIdx.x] = 0;

    /* Iterate over all FSRs */
    while (tid < *num_FSRs) {

	curr_material = &materials[FSR_materials[tid]];
	nu_sigma_f = curr_material->_nu_sigma_f;
	volume = FSR_volumes[tid];

	/* Iterate over all energy groups and update
	 * fission source for this block */
	for (int e=0; e < *num_groups; e++) {
	    source = nu_sigma_f[e] * scalar_flux(tid,e) * volume;
	    shared_fission_source[threadIdx.x] += source;
	}
		 
	/* Increment thread id */
	tid += blockDim.x * gridDim.x;
    }

    /* Copy this threads fission source to global memory */
    tid = threadIdx.x + blockIdx.x * blockDim.x;
    fission_sources[tid] = shared_fission_source[threadIdx.x];
    
    return;
}


/**
 * @brief Normalizes all flatsourceregion scalar fluxes and track boundary
 *        angular fluxes to the total fission source (times \f$ \nu \f$).
 * @param scalar_flux an array of the flat source region scalar fluxes
 * @param boundary_flux an array of the boundary fluxes
 * @param norm_factor the normalization factor
 */
__kernel void normalizeFluxesOnDevice(
        __global FP_PRECISION* scalar_flux, 
        __global FP_PRECISION* boundary_flux, 
        __global FP_PRECISION norm_factor) {

    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    
    /* Normalize scalar fluxes for each flat source region */
    while(tid < *num_FSRs) {

        for (int e=0; e < *num_groups; e++)
	    scalar_flux(tid,e) *= norm_factor;

	tid += blockDim.x * gridDim.x;
    }

    tid = threadIdx.x + blockIdx.x * blockDim.x;

    /* Normalize angular boundary fluxes for each track */
    while(tid < *tot_num_tracks) {
        for (int pe2=0; pe2 < 2*(*polar_times_groups); pe2++)
	    boundary_flux(tid,pe2) *= norm_factor;

	tid += blockDim.x * gridDim.x;
    }

    return;
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
 *                    \left(\frac{Q^i - Q^{i-1}{Q^i}\right)^2}{# FSRs}}} \f$
 *
 * @param FSR_materials an array of flat source region material UIDs
 * @param materials an array of material pointers
 * @param scalar_flux an array of flat source region scalar fluxes
 * @param source an array of flat source region sources
 * @param reduced_source an array of flat source region sources / total xs
 * @param inverse_k_eff the inverse of keff
 * @param an array of the source residuals 
 * @return the residual between this source and the previous source
 */
__kernel void computeFSRSourcesOnDevice(
        __global int* FSR_materials,
        __global dev_material* materials,
        __global FP_PRECISION* scalar_flux,
        __global FP_PRECISION* source,
        __global FP_PRECISION* old_source,
        __global FP_PRECISION* reduced_source,
        __global FP_PRECISION inverse_k_eff,
        __global FP_PRECISION* source_residuals) {

    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    /* Reset the residual for the old and new fission sources to zero */
    source_residuals[threadIdx.x + blockIdx.x * blockDim.x] = 0.0;

    FP_PRECISION fission_source;
    FP_PRECISION scatter_source;

    dev_material* curr_material;

    double* nu_sigma_f;
    double* sigma_s;
    double* sigma_t;
    double* chi;

    /* Iterate over all FSRs */
    while (tid < *num_FSRs) {

	curr_material = &materials[FSR_materials[tid]];

	nu_sigma_f = curr_material->_nu_sigma_f;
	sigma_s = curr_material->_sigma_s;
	sigma_t = curr_material->_sigma_t;
	chi = curr_material->_chi;

	/* Initialize the fission source to zero for this region */
	fission_source = 0;
	
	/* Compute total fission source for current region */
	for (int e=0; e < *num_groups; e++)
	    fission_source += scalar_flux(tid,e) * nu_sigma_f[e];
      
	/* Compute total scattering source for region for group G */
	for (int G=0; G < *num_groups; G++) {
	    scatter_source = 0;
	
	    for (int g=0; g < *num_groups; g++)
	        scatter_source += 
		    sigma_s[G*(*num_groups)+g] * scalar_flux(tid,g);
	
	    /* Set the total source for this region in this group */
	    source(tid,G) = (inverse_k_eff * fission_source * chi[G] +
			     scatter_source) * ONE_OVER_FOUR_PI;

	    reduced_source(tid,G) = __fdividef(source(tid,G), sigma_t[G]);

	    /* Compute the norm of residuals of the sources for convergence */
	    if (fabs(source(tid,G)) > 1E-10)
	        source_residuals[threadIdx.x + blockIdx.x * blockDim.x] +=
		    pow((source(tid,G) - old_source(tid,G)) /
		         source(tid,G), 2);

	    /* Update the old source */	
	    old_source(tid,G) = source(tid,G);
	}
	
	/* Increment the thread id */
	tid += blockDim.x * gridDim.x;
    }

    return;
}

/**
 * @brief Compute the total fission source from all flat source regions.
 * @param FSR_volumes an array of the flat source region volumes
 * @param FSR_materials an array of the flat source region material UIDs
 * @param materials an array of the material pointers
 * @param scalar_flux an array of flat source region scalar fluxes
 * @param tot_absorption an array of flat source region absorption rates
 * @param tot_fission an array of flat source region fission rates
 */
__kernel void computeFissionAndAbsorption(
        __global FP_PRECISION* FSR_volumes,
        __global int* FSR_materials,
        __global dev_material* materials,
        __global FP_PRECISION* scalar_flux,
        __global FP_PRECISION* tot_absorption,
        __global FP_PRECISION* tot_fission) {

    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    dev_material* curr_material;
    double* nu_sigma_f;
    double* sigma_a;
    FP_PRECISION volume;

    double absorption = 0.;
    double fission = 0.;

    /* Iterate over all FSRs */
    while (tid < *num_FSRs) {
        
	curr_material = &materials[FSR_materials[tid]];
	nu_sigma_f = curr_material->_nu_sigma_f;
	sigma_a = curr_material->_sigma_a;
	volume = FSR_volumes[tid];

	double curr_abs = 0.;
	double curr_fission = 0.;

	/* Iterate over all energy groups and update
	 * fission and absorption rates for this block */
	for (int e=0; e < *num_groups; e++) {
	    curr_abs += sigma_a[e] * scalar_flux(tid,e);
	    curr_fission += nu_sigma_f[e] * scalar_flux(tid,e);
	}

	absorption += curr_abs * volume;
	fission += curr_fission * volume;

	/* Increment thread id */
	tid += blockDim.x * gridDim.x;
    }

    /* Copy this thread's fission and absorption rates to global memory */
    tid = threadIdx.x + blockIdx.x * blockDim.x;
    tot_absorption[tid] = absorption;
    tot_fission[tid] = fission;

    return;
}


/**
 * @brief Compute the index into the exponential prefactor hashtable.
 * @details This method computes the index into the exponential prefactor
 *          hashtable for a segment length multiplied by the total 
 *          cross-section of the material the segment resides in.
 * @param sigm_t_l the cross-section multiplied by segment length
 * @return the hasthable index
 */ 

// TODO: is there separate syntax for helper functions?
__kernel int computePrefactorIndex(FP_PRECISION tau) {
    int index = tau * *inverse_prefactor_spacing;
    index *= *two_times_num_polar;
    return index;
}


/**
 * @brief Perform an atomic addition in double precision to an array address.
 * @details This method is straight out of CUDA C Developers Guide (cc 2013)
 * @param address the array memory address
 * @param value the value to add to the array 
 * @return the atomically added array value and input value
 *
 */

// TODO: is there separate syntax for helper functions?
__kernel double atomicAdd(
    double* address,
    double val) {

    unsigned long long int* address_as_ull = (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;

    do {
        assumed = old;
	old = atomicCAS(address_as_ull, assumed, __double_as_longlong(val +
			__longlong_as_double(assumed)));
    } while (assumed != old);
  
    return __longlong_as_double(old);
}


/**
 * @brief Computes the exponential term in the transport equation for a
 *        track segment.
 * @details This method computes \f$ 1 - exp(-l\Sigma^T_g/sin(\theta_p)) \f$ 
 *          for a segment with total group cross-section and for
 *          some polar angle.
 * @brief sigma_t the total group cross-section at this energy
 * @brief length the length of the line segment projected in the xy-plane
 * @param _prefactor_array the exponential prefactor interpolation table
 * @brief p the polar angle index
 * @return the evaluated exponential
 */
__device__ FP_PRECISION computeExponential(
        FP_PRECISION sigma_t, 
        FP_PRECISION length,
        __global FP_PRECISION* _prefactor_array,
        int p) {

    FP_PRECISION exponential;
    FP_PRECISION tau = sigma_t * length;

    /* Evaluate the exponential using the lookup table - linear interpolation */
    if (*interpolate_exponential) {
        int index;

	index = int(tau * (*inverse_prefactor_spacing)) * (*two_times_num_polar);
	exponential = (1. - (_prefactor_array[index+2 * p] * tau + 
			  _prefactor_array[index + 2 * p +1]));
    }

    /* Evalute the exponential using the intrinsic exp function */
    else {
        FP_PRECISION sintheta = sinthetas[p];
	#ifdef SINGLE
	exponential = 1.0 - __expf(- tau / sintheta);
	#else
	exponential = 1.0 - exp(- tau / sintheta);
	#endif
    }

    return exponential;
}


/**
 * @brief Computes the contribution to the flat source region scalar flux
 *        from a single track segment and a single energy group.
 * @details This method integrates the angular flux for a track segment across
 *        energy groups and polar angles, and tallies it into the flat
 *        source region scalar flux, and updates the track's angular flux.
 * @param curr_segment a pointer to the segment of interest
 * @param energy_group the energy group of interest
 * @param materials the array of materials
 * @param track_flux a pointer to the track's angular flux
 * @param reduced_source the array of flat source region sources / total xs
 * @param polar_weights the array of polar quadrature weights
 * @param _prefactor_array the exponential prefactor interpolation table
 * @param scalar_flux the array of flat source region scalar fluxes
 */
__kernel void scalarFluxTally(
        __global dev_segment* curr_segment, 
        int energy_group,
        __global dev_material* materials,
        __global FP_PRECISION* track_flux,
        __global FP_PRECISION* reduced_source,
        __global FP_PRECISION* polar_weights,
        __global FP_PRECISION* _prefactor_array,
        __global FP_PRECISION* scalar_flux) {

    int fsr_id = curr_segment->_region_uid;
    FP_PRECISION length = curr_segment->_length;
    dev_material* curr_material = &materials[curr_segment->_material_uid];
    double *sigma_t = curr_material->_sigma_t;

    /* The average flux long this segment in this flat source region */
    FP_PRECISION psibar;
    FP_PRECISION exponential;

    /* Zero the FSR scalar flux contribution from this segment 
     * and energy group */
    FP_PRECISION fsr_flux = 0.0;
    
    /* Compute the exponential prefactor hashtable index */
    
    /* Loop over polar angles */
    for (int p=0; p < *num_polar; p++) {
        exponential = computeExponential(sigma_t[energy_group], 
					 length, _prefactor_array, p);
        psibar = (track_flux[p] - reduced_source(fsr_id,energy_group)) * 
	         exponential;
	fsr_flux += psibar * polar_weights[p];
	track_flux[p] -= psibar;
    }

    /* Atomically increment the scalar flux for this flat source region */
    atomicAdd(&scalar_flux(fsr_id,energy_group), fsr_flux);
}


/**
 * @brief Updates the boundary flux for a track given boundary conditions.
 * @details For reflective boundary conditions, the outgoing boundary flux
 *          for the track is given to the reflecting track. For vacuum
 *          boundary conditions, the outgoing flux tallied as leakage.
 *          NOTE: Only one energy group is transferred by this routine.
 * @param curr_track a pointer to the track of interest
 * @param track_flux an array of the outgoing track flux
 * @param boundary_flux an array of all angular fluxes
 * @param leakage an array of leakages for each CUDA thread
 * @param polar_weights an array of polar quadrature weights
 * @param energy_angle_index the energy group index
 * @param direction the track direction (forward - true, reverse - false)
 */
__kernel void transferBoundaryFlux(
        __global dev_track* curr_track,
        __global FP_PRECISION* track_flux, 
        __global FP_PRECISION* boundary_flux, 
        __global FP_PRECISION* leakage,
        __global FP_PRECISION* polar_weights,
        __global int energy_angle_index,
        __global bool direction) {

    int start = energy_angle_index;
    bool bc;
    int track_out_id;

    /* Extract boundary conditions for this track and the pointer to the
     * outgoing reflective track, and index into the leakage array */

    /* For the "forward" direction */
    if (direction) {
        bc = curr_track->_bc_out;
	track_out_id = curr_track->_track_out;
	start += curr_track->_refl_out * (*polar_times_groups);
    }

    /* For the "reverse" direction */
    else {
        bc = curr_track->_bc_in;
	track_out_id = curr_track->_track_in;
        start += curr_track->_refl_in * (*polar_times_groups);
    }
    
    FP_PRECISION* track_out_flux = &boundary_flux(track_out_id,start);

    /* Put track's flux in the shared memory temporary flux array */
    for (int p=0; p < *num_polar; p++) {
	    track_out_flux[p] = track_flux[p] * bc;
	    leakage[0] += track_flux[p] * polar_weights[p] * (!bc);
    }
}


/**
 * @brief This method performs one transport sweep of one halfspace of all 
 *        azimuthal angles, tracks, segments, polar angles and energy groups.
 * @details The method integrates the flux along each track and updates the 
 *          boundary fluxes for the corresponding output track, while updating 
 *          the scalar flux in each flat source region
 * @param scalar_flux an array of flat source region scalar fluxes
 * @param boundary_flux an array of boundary fluxes
 * @param reduced_source an array of flat source region sources / total xs
 * @param leakage an array of angular flux leakaages
 * @param materials an array of material pointers
 * @param tracks an array of tracks
 * @param _prefactor_array an array for the exponential prefactor table
 * @param tid_offset the track offset for azimuthal angle halfspace
 * @param tid_max the upper bound on the track IDs for this azimuthal 
 *                angle halfspace
 */
__kernel void transportSweepOnDevice(
        __global FP_PRECISION* scalar_flux,
        __global FP_PRECISION* boundary_flux,
        __global FP_PRECISION* reduced_source,
        __global FP_PRECISION* leakage,
        __global dev_material* materials,
        __global dev_track* tracks,
        __global FP_PRECISION* _prefactor_array,
        int tid_offset,
        int tid_max) {

    /* Shared memory buffer for each thread's angular flux */
    extern __local FP_PRECISION temp_flux[];
    FP_PRECISION* track_flux;

    int tid = tid_offset + threadIdx.x + blockIdx.x * blockDim.x;

    int track_id = int(tid / *num_groups);
    int track_flux_index = threadIdx.x * (*two_times_num_polar);
    int energy_group = tid % (*num_groups);
    int energy_angle_index = energy_group * (*num_polar);

    dev_track* curr_track;
    int num_segments;
    dev_segment* curr_segment;

    /* Iterate over track with azimuthal angles in (0, pi/2) */
    while (track_id < tid_max) {

        /* Initialize local registers with important data */
        curr_track = &tracks[track_id];
      	num_segments = curr_track->_num_segments;

	/* Retrieve a pointer to this thread's shared memory buffer for angular flux */
	track_flux = &temp_flux[track_flux_index];
      
	/* Put track's flux in the shared memory temporary flux array */
      	for (int p=0; p < *num_polar; p++) {
	
	    /* Forward flux along this track */
	    track_flux[p] = boundary_flux(track_id,p+energy_angle_index);

	    /* Reverse flux along this track */
	    track_flux[(*num_polar) + p] = boundary_flux(track_id,p+energy_angle_index+(*polar_times_groups));
      	}
      
	/* Loop over each segment in forward direction */
	for (int i=0; i < num_segments; i++) {
	    curr_segment = &curr_track->_segments[i];
	    scalarFluxTally(curr_segment, energy_group, materials,
			    track_flux, reduced_source, polar_weights,
			    _prefactor_array, scalar_flux);
	}

	/* Transfer flux to outgoing track */
	transferBoundaryFlux(curr_track, track_flux, boundary_flux, 
			     &leakage[threadIdx.x + blockIdx.x * blockDim.x], 
			     polar_weights, energy_angle_index, true);

	/* Loop over each segment in reverse direction */
	track_flux = &temp_flux[track_flux_index + (*num_polar)];

	for (int i=num_segments-1; i > -1; i--) {
	    curr_segment = &curr_track->_segments[i];
	    scalarFluxTally(curr_segment, energy_group, materials,
			    track_flux, reduced_source, polar_weights,
			    _prefactor_array, scalar_flux);
	}

	/* Transfer flux to outgoing track */
	transferBoundaryFlux(curr_track, track_flux, boundary_flux, 
			     &leakage[threadIdx.x + blockIdx.x * blockDim.x], 
			     polar_weights, energy_angle_index, false);

        /* Update the indices for this thread to the next track, energy group */
	tid += blockDim.x * gridDim.x;
        track_id = int(tid / *num_groups);
	energy_group = tid % (*num_groups);
	energy_angle_index = energy_group * (*num_polar);
    }

    return;
}


/**
 * @brief Add the source term contribution in the transport equation to 
 *        the flat source region scalar flux
 * @param scalar_flux an array of flat source region scalar fluxes
 * @param reduced_source an array of flat source region sources / total xs
 * @param FSR_volumes an array of flat source region volumes
 * @param FSR_materials an array of flat source region material UIDs
 * @param materials an array of material pointers
 */
__kernel void addSourceToScalarFluxOnDevice(
        __global FP_PRECISION* scalar_flux,
        __global FP_PRECISION* reduced_source,
        __global FP_PRECISION* FSR_volumes,
        __global int* FSR_materials,
        __global dev_material* materials) {
  
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    FP_PRECISION volume;
    
    dev_material* curr_material;
    double* sigma_t;

    /* Iterate over all FSRs */
    while (tid < *num_FSRs) {

	curr_material = &materials[FSR_materials[tid]];
	volume = FSR_volumes[tid];
	sigma_t = curr_material->_sigma_t;
	
	/* Iterate over all energy groups */
	for (int i=0; i < *num_groups; i++) {
	    scalar_flux(tid,i) *= 0.5;
	    scalar_flux(tid,i) = FOUR_PI * reduced_source(tid,i) + 
	      __fdividef(scalar_flux(tid,i), (sigma_t[i] * volume));
	}

	/* Increment thread id */
	tid += blockDim.x * gridDim.x;
    }

    return;
}


/**
 * @brief Computes the volume-weighted, energy integrated fission rate in 
 *        each flat source region and stores them in an array indexed by 
 *        flat source region ID.
 * @param fission_rates an array to store the fission rates, passed in as a
 *        numpy array from Python
 * @param fission_rates an array in which to store the flat source region fission rates
 * @param FSR_materials an array of flat source region material UIDs
 * @param materials an array of material pointers
 * @param scalar_flux an array of flat source region scalar fluxes
 */
__kernel void computeFSRFissionRatesOnDevice(
        __global double* fission_rates,
        __global int* FSR_materials,
        __global dev_material* materials,
        __global FP_PRECISION* scalar_flux) {

    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    dev_material* curr_material;
    double* sigma_f;

    /* Loop over all FSRs and compute the volume-weighted fission rate */
    while (tid < *num_FSRs) {

	curr_material = &materials[FSR_materials[tid]];
	sigma_f = curr_material->_sigma_f;

	/* Initialize the fission rate for this FSR to zero */
	fission_rates[tid] = 0.0;

        for (int i=0; i < *num_groups; i++)
	    fission_rates[tid] += sigma_f[i] * scalar_flux(tid,i);

	/* Increment thread id */
	tid += blockDim.x * gridDim.x;
    }

    return;
}
