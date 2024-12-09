#ifndef ROTATION_H
#define ROTATION_H

#include "gsl_quaternion_float.h"
#include "gsl/gsl_vector_float.h"
#include "gsl/gsl_matrix_float.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_blas_types.h"


typedef gsl_matrix_float rotMat; 

typedef enum rot_method_e{
  QUAT_METHOD,
  ROTMAT_METHOD
} rot_method_t;

typedef struct rotation_s {
  union rot_u {
    rotMat *mat;
    gsl_quat_float *quat_invPair[2]; //quat_invPair[0] *q | quat_invPair[1] *(q^-1) 
  } value;
  rot_method_t method;
} rotation_t;

extern int createRotationFromQuat(gsl_quat_float *pQ, rotation_t *rot);

extern int rotateVector(gsl_vector_float *v, rotation_t *rot);

#endif