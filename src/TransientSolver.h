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
#include "CmfdTransient.h"
#include "Mesh.h"
#include "Solver.h"
#include "Timer.h"
#include "TimeStepper.h"
#ifdef _mm_malloc
#undef _mm_malloc
#undef _mm_free
#endif
#include "petsc.h"
#include <petscmat.h>


class TransientSolver {
private:
  Geometry* _geom;
  CmfdTransient* _cmfd;
  Solver* _solver;
  TimeStepper* _ts;
  Material** _materials;
  Mesh* _mesh;
  Timer* _timer;
  
  bool _implicit;
  bool _improved;
  bool _steady_state;
  bool _adj_weight;
  bool _multigroup_PKE;
  double _dt_outer;
  double _dt_intermediate;
  double _dt_pke;
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
  std::vector<double> _step_times;
  std::vector<double> _dts_outer;
  std::vector<double> _dts_middle;
  int _num_pke_groups;
  int _num_delay_groups;
  int _num_groups;
  
  double _nu;
  double _kappa;
  double _alpha;
  
  double* _N;
  double* _B;
  double* _B_prev;
  double* _B_inter;
  double* _C;
  
  
public:
  TransientSolver(Geometry* geom, CmfdTransient* cmfd,
		  Solver* solver=NULL);
  virtual ~TransientSolver();
  
  double computePower(materialState state=CURRENT);
  void updateTemperatures();
  double computeResidual();
  void printProblemSpecs();
  double computeNewShape(materialState state=CURRENT);
  void trimVectors(int len_conv);
  void initializePrecursorConc();
  void updatePrecursorConc(materialState state);
  void copyAmplitude(double* b_from, double* b_to);
  void setStartTime(double time);
  void setEndTime(double time);
  void initializeTimeStepper();
  void initializeTransientMaterials();
  void setAdjointWeighting(bool weighting);
  void normalizeFlux(materialState state=CURRENT);
  
  void solveInitialState();
  void solveOuterStep();
  
  void constructN();
  
  void setDtOuter(double dt);
  void setDtIntermediate(double dt);
  void setDtPKE(double dt);
  
  void setKappa(double kappa);
  void setNu(double nu);
  void setAlpha(double alpha);
  void setNumDelayGroups(int num_groups);
  
  void updateFlux();
  double computeCoreTemp();
  
  void computeVolCore();
  void setTransientMethod(const char* trans_type);
  void setOuterImplicit(bool implicit);
  void setImproved(bool improved);
  void setMultigroupPKE(bool multigroup);
  void setSteadyState(bool steady_state);
  void setPowerInit(double power);
  void syncMaterials(materialState state);
  void copyPrecConc(materialState state_from, materialState state_to);
  void copyTemperature(materialState state_from, materialState state_to);
  void copyFieldVariables(materialState state_from, materialState state_to);
};


#endif /* TRANSIENTSOLVER_H_ */
