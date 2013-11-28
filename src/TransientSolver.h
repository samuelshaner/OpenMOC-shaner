/*
 * TransientSolver.h
 *
 *  Created on: Dec 15, 2012
 *      Author: samuelshaner
 */

#ifndef TRANSIENTSOLVER_H_
#define TRANSIENTSOLVER_H_

#include <math.h>
#include <string>
#include <sstream>
#include <vector>
#include <ctime>
#include <time.h>
#include "Geometry.h"
#include "Quadrature.h"
#include "log.h"
#include "Tcmfd.h"
#include "Mesh.h"
#include "Cmfd.h"
#include "ThreadPrivateSolverTransient.h"
#include "Timer.h"
#include "TimeStepper.h"

class TransientSolver {
private:
    Geometry* _geom;
    Tcmfd* _tcmfd;
    Cmfd* _cmfd;
    Solver* _solver;
    TimeStepper* _ts;
    Mesh* _mesh;
    Mesh* _geom_mesh;
    Timer* _timer;
    
    double _dt_moc;
    double _dt_cmfd;
    double _power_init;
    double _vol_core;
    double _k_eff_0;
    double _start_time;
    double _end_time;
    std::vector<double> _temp;
    std::vector<double> _power;
    std::vector<double> _time;
    double _power_factor;
    solveType _solve_method;
    transientType _transient_method;
    int _ndg;
    int _ng;
    
    double _nu;
    double _kappa;
    double _alpha;
    
    
public:
    TransientSolver(Geometry* geom, Tcmfd* tcmfd, Cmfd* cmfd,
		    Solver* solver=NULL);
    virtual ~TransientSolver();
    
    /* worker functions */
    double computeCoreTemp();  
    void computeVolCore();
    void sync(materialState state);
    void copyPrecConc(materialState state_from, materialState state_to);
    void copyTemperature(materialState state_from, materialState state_to);
    void copyFieldVariables(materialState state_from, materialState state_to);
    double computePower(materialState state=CURRENT);
    void updateTemperatures();
    double computeResidual(solveType solve_type);
    void trimVectors(int len_conv);
    void initializePrecursorConc();
    void updatePrecursorConc();
    void initializeTimeStepper();
    void initializeTransientMaterials();
    void solveInitialState();
    void solveOuterStep();
    void mapPrecConc();
    
    /* setters */
    void setDtMOC(double dt);
    void setDtCMFD(double dt);  
    void setKappa(double kappa);
    void setNu(double nu);
    void setAlpha(double alpha);
    void setNumDelayGroups(int num_groups);
    void setTransientMethod(const char* trans_type);
    void setPowerInit(double power);
    void setStartTime(double time);
    void setEndTime(double time);
    double getPower();
    double getTime();
    double getTemp();



};

#endif /* TRANSIENTSOLVER_H_ */
