#ifndef JACOBIAN_H
#define JACOBIAN_H

#include "gsl/gsl_vector_float.h"
#include "gsl/gsl_matrix_float.h"
#include "inttypes.h"

typedef void (*diff_function)(gsl_vector_float* in, gsl_vector_float* out);

extern void jacobian(gsl_matrix_float *pJacMatrix, diff_function pFunc, gsl_vector_float *pVec0, float delta);

#endif