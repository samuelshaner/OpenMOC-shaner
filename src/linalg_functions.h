/*
 * linalg_functions.h
 *
 *  Created on: November 11, 2013
 *  Author: Samuel Shaner
 */

#ifndef LINALGFUNCTIONS_H_
#define LINALGFUNCTIONS_H_

#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include <omp.h>

/* function definitions */
inline double rayleighQuotient(double* A, double* M, double* x, double* snew, double* sold, int cx, int cy, int ng); 
inline void matSubtract(double* AM, double* A, double omega, double* M, int cx, int cy, int ng);
inline double vecMax(double* vec, int length);
inline void vecScale(double* vec, double scale_val, int length);
inline void vecSet(double* vec, double val, int length);
inline void matMultM(double* mat, double* vec_x, double* vec_y, int num_cells, int num_groups);
inline void matMultA(double* mat, double* vec_x, double* vec_y, int cells_x, int cells_y, int num_groups);
inline void vecSubtract(double* a, double* b, double* c, int ng);
inline double vecSum(double* vec, int length);
inline void vecCopy(double* vec_from, double* vec_to, int length);
inline void vecZero(double* vec, int length);
inline double norm(double* x, int nc);
inline void vecWAXPY(double* vec_w, double a, double* vec_x, double* vec_y, int length);
inline void vecWAXPY(float* vec_w, double a, double* vec_x, double* vec_y, int length);
inline void linearSolve(double* mat, double* vec_x, double* vec_b, double* vec_x_old, double conv, double omega, int cx, int cy, int ng, int max_iter);
inline void dumpVector(double* vec, int length);

double rayleighQuotient(double* A, double* M, double* x, double* snew, double* sold, int cx, int cy, int ng){

    double numer = 0.0;
    double denom = 0.0;
    
    matMultA(A, x, sold, cx, cy, ng);
    matMultM(M, x, snew, cx*cy, ng);
    
    for (int i = 0; i < cx*cy*ng; i++){
	numer += x[i]*snew[i];
	denom += x[i]*sold[i];
    }

    return numer/denom;
}


void matSubtract(double* AM, double* A, double omega, double* M, int cx, int cy, int ng){

    /* copy A to AM */
    vecCopy(A, AM, cx*cy*ng*(ng+4));
    
    /* subtract omega*M from AM */
    for (int i = 0; i < cx*cy*ng; i++){
	for (int e = 0; e < ng; e++)
	    AM[i*(ng+4)+e+2] -= omega*M[i*ng+e];		
    }
}


double vecMax(double* vec, int length){
    
    double max = vec[0];
    
    for (int i = 0; i < length; i++)
	max = std::max(max, vec[i]);

    return max;
}


void vecScale(double* vec, double scale_val, int length){

    for (int i = 0; i < length; i++)
      vec[i] *= scale_val;
}


void vecSet(double* vec, double val, int length){

    for (int i = 0; i < length; i++)
	vec[i] = val;
}


void vecSubtract(double* a, double* b, double* c, int ng){

    for (int i = 0; i < ng; i++)
	c[i] = a[i] - b[i];

}


double norm(double* x, int nc){

    double norm = 0.0;
    
    for (int i = 0; i < nc; i++)
	norm += pow(x[i], 2);

    norm = pow(norm, 0.5);
    return norm;
}


void matMultM(double* mat, double* vec_x, double* vec_y, int num_cells, int num_groups){

    vecZero(vec_y, num_cells*num_groups);

    for (int i = 0; i < num_cells; i++){
	for (int g = 0; g < num_groups; g++){
	    for (int e = 0; e < num_groups; e++){
		vec_y[i*num_groups+g] += mat[(i*num_groups+g)*num_groups+e] * vec_x[i*num_groups+e];
	    }
	}
    }
}


void matMultA(double* mat, double* vec_x, double* vec_y, int cells_x, int cells_y, int num_groups){

    vecZero(vec_y, cells_x*cells_y*num_groups);
    int row;
    
    for (int y = 0; y < cells_y; y++){
	for (int x = 0; x < cells_x; x++){
	    for (int g = 0; g < num_groups; g++){
		row = (y*cells_x+x)*num_groups + g;
		
		/* left surface */
		if (x != 0)
		    vec_y[row] += mat[row*(num_groups+4)] * vec_x[(y*cells_x+x-1)*num_groups +g]; 
		
		/* bottom surface */
		if (y != cells_y - 1)
		    vec_y[row] += mat[row*(num_groups+4)+1] * vec_x[((y+1)*cells_x+x)*num_groups +g]; 
		
		/* right surface */
		if (x != cells_x-1)
		    vec_y[row] += mat[row*(num_groups+4)+num_groups+2] * vec_x[(y*cells_x+x+1)*num_groups +g]; 
		
		/* top surface */
		if (y != 0)
		    vec_y[row] += mat[row*(num_groups+4)+num_groups+3] * vec_x[((y-1)*cells_x+x)*num_groups +g]; 
		
		for (int e = 0; e < num_groups; e++){
		    vec_y[row] += mat[row*(num_groups+4)+2+e] * vec_x[(y*cells_x+x)*num_groups+e];
		}
	    }
	}
    } 
}


double vecSum(double* vec, int length){
   
    double sum = 0.0;

    for (int i = 0; i < length; i++)
	sum += vec[i];

    return sum;
}


void vecCopy(double* vec_from, double* vec_to, int length){

  for (int i = 0; i < length; i++)
    vec_to[i] = vec_from[i];
}


void vecZero(double* vec, int length){

  for (int i = 0; i < length; i++)
    vec[i] = 0.0;
}


void vecWAXPY(double* vec_w, double a, double* vec_x, double* vec_y, int length){
    
    for (int i = 0; i < length; i++){
	vec_w[i] = a * vec_x[i] + vec_y[i];
    }    
}


void vecWAXPY(float* vec_w, double a, double* vec_x, double* vec_y, int length){
    
    for (int i = 0; i < length; i++){
	vec_w[i] = a * vec_x[i] + vec_y[i];
    }    
}


void linearSolve(double* mat, double* vec_x, double* vec_b, double* vec_x_old, double conv, double omega, int cx, int cy, int ng, int max_iter){

    double norm = 1e10;
    int row = 0;
    double val = 0.0;
    int iter = 0;

    /* perform GS iteration */
    while (norm > conv){

	/* pass new flux to old flux */
	vecCopy(vec_x, vec_x_old, cx*cy*ng);

	for (int y = 0; y < cy; y++){
	    for (int x = 0; x < cx; x++){
		for (int g = 0; g < ng; g++){
		    row = ((y*cx+x)*ng+g);
		    val = 0.0;
		    
		    val += (1.0 - omega) * vec_x[row];
		    
		    /* source */
		    val += omega * vec_b[row] / mat[row*(ng+4)+g+2];
		    
		    /* left surface */
		    if (x != 0)
			val -= omega * vec_x[(y*cx+x-1)*ng+g] * mat[row*(ng+4)] / mat[row*(ng+4)+g+2];
		   	    
		    /* bottom surface */
		    if (y != cy - 1)
			val -= omega * vec_x[((y+1)*cx+x)*ng+g] * mat[row*(ng+4)+1] / mat[row*(ng+4)+g+2];
		    
		    /* group to group */
		    for (int e = 0; e < ng; e++){
			if (e != g)
			    val -= omega * vec_x[(y*cx+x)*ng+e] * mat[row*(ng+4)+2+e] / mat[row*(ng+4)+2+g];
		    }
		    
		    /* right surface */
		    if (x != cx - 1)
			val -= omega * vec_x[(y*cx+x+1)*ng+g] * mat[row*(ng+4)+ng+2] / mat[row*(ng+4)+g+2];
		    
		    /* top surface */
		    if (y != 0)
			val -= omega * vec_x[((y-1)*cx+x)*ng+g] * mat[row*(ng+4)+ng+3] / mat[row*(ng+4)+g+2];
		    
		    vec_x[row] = val;
		}
	    }
	}
	
	norm = 0.0;
	for (int i = 0; i < cx*cy*ng; i++){
	    if (vec_x[i] != 0.0)
		norm += pow((vec_x[i] - vec_x_old[i])/(vec_x[i]), 2);
	}

	norm = pow(norm, 0.5) / (cx*cy*ng);
	iter++;

        //std::cout << "vec x sum: " << vecSum(vec_x, cx*cy*ng) << std::endl;

	if (iter == max_iter)
	    break;
    }

    //std::cout << "GS iterations: " << iter << std::endl;
}


void dumpVector(double* vec, int length){

    for (int i = 0; i < length; i++)
	std::cout << "vec i: " << i << ", value: " << vec[i] << std::endl;

}


#endif /* LINALG_FUNCTIONS_H_ */
