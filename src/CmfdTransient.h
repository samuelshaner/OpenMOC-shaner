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

  double* _frequency;
  double* _amplitude;
  transientType _transient_method;
  double _beta_sum;
  int _num_delay_groups;
  bool _multigroup_PKE;
  bool _adj_weight;
  bool _initial_solve;
  
  /* materials parameters */
  double* _beta;
  double* _lambda;
  double* _velocity;
  
public:
  CmfdTransient(Geometry* geometry, double criteria=1e-8);
  virtual ~CmfdTransient();
  virtual void constructMatrices();
  virtual double computeKeff(fluxType flux_method=PRIMAL);
  void setTransientType(transientType trans_type);
  void setFrequency(int i, double value);
  void setAmplitude(int i, double value);
  void setKeff0(double keff_0);
  virtual void rescaleFlux();
  void setMultigroupPKE(bool mg);
  transientType getTransientType();
  void initialize(TimeStepper* _time_stepper, bool adjoint_weighted);
  
  void setBeta(double* beta, int num_delay_groups);
  void setLambda(double* decay_const, int num_delay_groups);
  void setVelocity(double* velocity, int num_groups);
  
  double* getBeta();
  double* getLambda();
  double* getVelocity();
  void setInitialSolve(bool initial);
  
  void matMultAdd(double* mat, double* vec_1, double* vec_2, double* vec_3);
  void vecWAXPY(double* vec_w, double a, double* vec_x, double* vec_y);

};

#endif /* CMFD_TRANSIENT_H_ */
