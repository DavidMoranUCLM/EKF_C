#ifndef MATH_UTILS_H
#define MATH_UTILS_H

#include "gsl/gsl_vector_float.h"
#include "gsl/gsl_matrix_float.h"

/**
 * @brief Creates a skew symetric matrix from a vector
 * 
 * @param pV Pointer of the source vector
 * @param pM Pointer to the result matrix
 */
extern int skewSymFromVector(gsl_vector_float *pV, gsl_matrix_float *pM);

extern void cross_product(const gsl_vector_float *a,
                               const gsl_vector_float *b, gsl_vector_float *c);
#endif