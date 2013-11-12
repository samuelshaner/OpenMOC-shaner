/**
 * @file Tcmfd.h
 * @brief The Tcmfd class.
 * @date October 14, 2013
 * @author Sam Shaner, MIT, Course 22 (shaner@mit.edu)
 */

#ifndef TCMFD_H_
#define TCMFD_H_

#ifdef __cplusplus
#define _USE_MATH_DEFINES
#include <math.h>
#include "Material.h"
#include "TimeStepper.h"
#include "FunctionalMaterial.h"
#include "Geometry.h"
#include "Timer.h"
#include "Mesh.h"
#include "linalg_functions.h"
#include <vector>
#endif

enum transientType {
  NONE,
  ADIABATIC,
  IQS,
  OMEGA_MODE
};


class Tcmfd{
protected:
  TimeStepper* _time_stepper;
  Mesh* _mesh;
  Geometry* _geom;
  Timer* _timer;

  /* integer values */
  int _cells_x;
  int _cells_y;
  int _num_groups;
  int _nc;
  int _num_delay_groups;

  /* float values */
  double _omega;
  double _k_eff_0;
  double _beta_sum;
  double _dt_cmfd;
  double _conv_criteria;
  
  /* matrix and vector objects */
  double* _b_prime;
  double* _b;
  double* _A;
  double* _M;
  double* _phi_old;
  double* _phi_new;
  double* _phi_temp;
  double* _snew;
  double* _sold;
 
  /* materials parameters */
  double* _beta;
  double* _lambda;
  double* _velocity;

  /* solve method (DIFFUSION or MOC) */
  solveType _solve_method; 
  transientType _transient_method;
  bool _initial_state; 
 
  
public:
  Tcmfd(Geometry* geometry, double criteria=1e-8);
  virtual ~Tcmfd();
  void constructMatrices();
  void solveTCMFD();
  void setTransientType(transientType trans_type);
  void setKeff0(double keff_0);
  transientType getTransientType();
  void setTimeStepper(TimeStepper* _time_stepper);
  
  void setBeta(double* beta, int num_delay_groups);
  void setLambda(double* decay_const, int num_delay_groups);
  void setVelocity(double* velocity, int num_groups);
  void setDtCMFD(double dt);

  double* getBeta();
  double* getLambda();
  double* getVelocity();
  
  void setInitialState(bool state);
  void setNumDelayGroups(int num_groups);

  Mesh* getMesh();
  void setOmega(double omega);
};

#endif /* TCMFD_H_ */
