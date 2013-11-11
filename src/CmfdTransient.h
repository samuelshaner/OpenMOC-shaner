/**
 * @file Cmfd.h
 * @brief The Cmfd class.
 * @date October 14, 2013
 * @author Sam Shaner, MIT, Course 22 (shaner@mit.edu)
 */

#ifndef CMFDTRANSIENT_H_
#define CMFDTRANSIENT_H_

#ifdef __cplusplus
#define _USE_MATH_DEFINES
#include <math.h>
#include "Material.h"
#include "TimeStepper.h"
#include "FunctionalMaterial.h"
#include "Cmfd.h"
#include <vector>
#endif

enum transientType {
  NONE,
  ADIABATIC,
  IQS,
  OMEGA_MODE
};


class CmfdTransient : public Cmfd {
protected:
  TimeStepper* _time_stepper;
  double* _b_prime;
  double* _b;

  double _k_eff_0;

  transientType _transient_method;
  double _beta_sum;
  int _num_delay_groups;
  
  /* materials parameters */
  double* _beta;
  double* _lambda;
  double* _velocity;

  double _dt_cmfd;
  bool _initial_state;

  
public:
  CmfdTransient(Geometry* geometry, double criteria=1e-8);
  virtual ~CmfdTransient();
  virtual void constructMatrices();
  virtual double computeKeff(fluxType flux_method=PRIMAL);
  void setTransientType(transientType trans_type);
  void setKeff0(double keff_0);
  transientType getTransientType();
  void initialize(TimeStepper* _time_stepper);
  
  void setBeta(double* beta, int num_delay_groups);
  void setLambda(double* decay_const, int num_delay_groups);
  void setVelocity(double* velocity, int num_groups);
  void setDtCMFD(double dt);

  double* getBeta();
  double* getLambda();
  double* getVelocity();
  
  void vecWAXPY(double* vec_w, double a, double* vec_x, double* vec_y);

  void setInitialState(bool state);
  void setNumDelayGroups(int num_groups);
};

#endif /* CMFD_TRANSIENT_H_ */
