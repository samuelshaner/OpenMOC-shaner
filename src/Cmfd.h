/**
 * @file Cmfd.h
 * @brief The Cmfd class.
 * @date October 14, 2013
 * @author Sam Shaner, MIT, Course 22 (shaner@mit.edu)
 */

#ifndef CMFD_H_
#define CMFD_H_

#ifdef __cplusplus
#define _USE_MATH_DEFINES
#include <utility>
#include <math.h>
#include <unordered_map>
#include <limits.h>
#include <string>
#include <sstream>
#include <queue>
#include <iostream>
#include <fstream>
#include "Quadrature.h"
#include "log.h"
#include "Mesh.h"
#include "Material.h"
#include "Surface.h"
#include "Geometry.h"
#include "Timer.h"
#include "linalg_functions.h"
#endif

class Cmfd {
protected:

  Geometry* _geometry;
  Mesh* _mesh;
  Timer* _timer;  

  /* keff */
  double _k_eff;

  /* matrix and vector objects */
  double* _A;
  double* _M;
  double* _phi_old;
  double* _phi_new;
  double* _sold;
  double* _snew;
  double* _phi_temp;

  /* float values */
  double _conv_criteria;
  double _omega;

  /* integer values */
  int _cells_x;
  int _cells_y;
  int _num_groups;
  int _nc;

  /* solve method (DIFFUSION or MOC) */
  solveType _solve_method;
  
public:
	
  Cmfd(Geometry* geometry, double criteria=1e-8);
  virtual ~Cmfd();

  virtual void constructMatrices();
  double computeKeff();
  void rescaleFlux();
  double* getA();
  double* getM();
  Mesh* getMesh();
  double getKeff();
  void setOmega(double omega);
};

#endif /* CMFD_H_ */
