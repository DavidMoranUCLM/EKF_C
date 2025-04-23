#ifndef MAGCORRECTION_H
#define MAGCORRECTION_H


#include "gsl/gsl_blas.h"
#include "gsl/gsl_matrix_float.h"
#include "gsl/gsl_vector_float.h"

extern void correctMag(gsl_matrix_float* P1, gsl_vector_float* x1,
  const gsl_vector_float* mag, const gsl_vector_float *mag_sigma);
  

#endif /* MAGCORRECTION_H */