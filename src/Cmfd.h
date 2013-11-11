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
#endif

/**
 * Flux types
 */
enum fluxType {
	PRIMAL,
	ADJOINT
};

/** Indexing macro for the ghosted flux vectors */
#define _phi_ghost_new(x,y,e) (_phi_ghost_new[((y+1)*(_cells_x+2) + x + 1) * _num_groups + e])
#define _phi_ghost_old(x,y,e) (_phi_ghost_old[((y+1)*(_cells_x+2) + x + 1) * _num_groups + e])

class Cmfd {

protected:

  /* pointer to quadrature object */
  Quadrature* _quad;

  /* pointer to geometry object */
  Geometry* _geometry;

  /* pointer to mesh object */
  Mesh* _mesh;

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
  double* _phi_ghost_new;
  double* _phi_ghost_old;

  /* l2 norm */
  double _l2_norm;

  /* keff convergence criteria */
  double _conv_criteria;
  double _omega;

  /* num cells in x and y direction */
  int _cells_x;
  int _cells_y;

  /* pointer to timer object */
  Timer* _timer;

  /* flux type (PRIMAL or ADJOINT) */
  fluxType _flux_method;

  /* number of groups */
  int _num_groups;

  /* number of fsrs */
  int _num_fsrs;

  /* solve method (DIFFUSION or MOC) */
  solveType _solve_method;
  
  /* arrays for fsr parameters */
  FP_PRECISION* _FSR_volumes;
  Material** _FSR_materials;
  FP_PRECISION* _FSR_fluxes;
  
public:
	
  Cmfd(Geometry* geometry, double criteria=1e-8);
  virtual ~Cmfd();

  /* worker functions */
  virtual void constructMatrices();
  void computeDs();
  void computeXS();
  void updateMOCFlux();
  double computeDiffCorrect(double d, double h);
  virtual double computeKeff();
  void initializeFSRs();
  virtual void rescaleFlux();

  /* get arrays, values, and objects */
  double* getA();
  double* getM();
  Mesh* getMesh();
  double getKeff();

  /* set parameters */
  void setMeshCellFlux();

  /* set fsr parameters */
  void setFSRMaterials(Material** FSR_materials);
  void setFSRVolumes(FP_PRECISION* FSR_volumes);
  void setFSRFluxes(FP_PRECISION* scalar_flux);
  void setFluxType(const char* flux_type);

  void vecZero(double* my_vec);
  void matMultM(double* mat, double* vec_1, double* vec_2);
  void matMultA(double* mat, double* vec_1, double* vec_2);
  void linearSolve(double* mat, double* vec_x, double* vec_b, double conv, double omega);
  void vecScale(double* vec, double scale_val);
  double vecSum(double* vec);
  void vecCopy(double* vec_from, double* vec_to, int x_len, int y_len);
  void matZero(double* mat, int width);
  void dumpVec(double* vec, int length);
  void setOmega(double omega);
  void vecSet(double* vec, double val);
  void vecNormal(double* mat, double* vec);
};

#endif /* CMFD_H_ */
