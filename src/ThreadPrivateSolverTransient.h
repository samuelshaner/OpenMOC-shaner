/**
 * @file ThreadPrivateSolverTransient.h
 * @brief The ThreadPrivateSolverTransient class.
 * @date November 11, 2013
 * @author Samuel Shaner, MIT, Course 22 (shaner@mit.edu)
 */


#ifndef THREADPRIVATESOLVERTRANSIENT_H_
#define THREADPRIVATESOLVERTRANSIENT_H_

#ifdef __cplusplus
#define _USE_MATH_DEFINES
#include <math.h>
#include <omp.h>
#include <stdlib.h>
#include "CPUSolver.h"
#include "TimeStepper.h"
#endif

/** Indexing scheme for the thread private scalar flux for each thread in
 *  each flat source region and in each energy group */
#define _thread_flux(tid,r,e) (_thread_flux[(tid)*_num_FSRs*_num_groups+(r)*_num_groups+(e)])

/** Indexing scheme for the thread private surface currents for each thread in
 *  each flat source region and in each energy group */
#define _thread_currents(tid,r,e) (_thread_currents[(tid)*_num_mesh_cells*8*_num_groups+(r % _geometry->getMesh()->getNumCurrents())*_num_groups+(e)])

/**
 * @class ThreadPrivateSolverTransient ThreadPrivateSolverTransient.h "openmoc/src/ThreadPrivateSolverTransient.h"
 * @brief This is a subclass of the CPUSolver which uses thread private 
 *        arrays for the flat source region scalar fluxes to minimize OMP atomics.
 * @details Since this class stores a separate copy of the flat source region scalar
 *          fluxes for each OMP thread, the memory requirements are greater than for
 *          the CPUSolver.
 */
class ThreadPrivateSolverTransient : public CPUSolver {

protected:

    /** An array for the flat source region scalar fluxes for each thread */
    FP_PRECISION* _thread_flux;

    double* _thread_currents;
    transientType _transient_method;

    void initializeFluxArrays();

    void flattenFSRFluxes(FP_PRECISION value);
    void zeroSurfaceCurrents();
    void scalarFluxTally(segment* curr_segment, int azim_index, 
			 FP_PRECISION* track_flux,
			 FP_PRECISION* fsr_flux,
			 bool fwd);
    void scalarFluxPass(segment* curr_segment, int azim_index,
			 FP_PRECISION* track_flux,
			 FP_PRECISION* fsr_flux,
			 bool fwd);
    void reduceThreadScalarFluxes();
    void reduceThreadSurfaceCurrents();
    void transportSweep();
    void transportPass();
    virtual FP_PRECISION computeFSRSources();

public:
    ThreadPrivateSolverTransient(Geometry* geometry=NULL, 
			TrackGenerator* track_generator=NULL,
			Cmfd* cmfd=NULL);
    virtual ~ThreadPrivateSolverTransient();


    virtual FP_PRECISION convergeSource(int max_iterations);

    void initialize();
    void scaleTrackFlux(double scale_val);
    void resetSegmentMaterials();
};


#endif /* THREADPRIVATESOLVERTRANSIENT_H_ */
