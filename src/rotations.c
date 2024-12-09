#include "rotations.h"

/**
 *
 * Private definitions
 *
 */

/**
 *
 * Private function declarations
 *
 */

/**
 *
 * Public function definitions
 *
 */

int rotateVector(gsl_vector_float *v, rotation_t *rot) {
  switch (rot->method) {
    case QUAT_METHOD:

      gsl_quat_float *pQv = gsl_quat_float_alloc();
      gsl_quat_float *pQMid = gsl_quat_float_alloc();
      gsl_quat_float_fromVector(v, pQv);

      gsl_vector_float_memcpy(pQMid, rot->value.quat_invPair[0]);
      gsl_quat_float_product(pQMid, pQv);

      gsl_quat_float_product(pQMid, rot->value.quat_invPair[1]);

      
      gsl_vector_float_set(v,0,gsl_quat_float_get(pQMid, 1));
      gsl_vector_float_set(v,1,gsl_quat_float_get(pQMid, 2));
      gsl_vector_float_set(v,2,gsl_quat_float_get(pQMid, 3));

      gsl_quat_float_free(pQv);
      gsl_quat_float_free(pQMid);

      break;

    case ROTMAT_METHOD:
      return -1;
      break;

    default:
      return -1;
  }

  return 0;
}

int createRotationFromQuat(gsl_quat_float *pQ, rotation_t *rot){
  rot->method = QUAT_METHOD;

  rot->value.quat_invPair[0] = gsl_quat_float_alloc();
  rot->value.quat_invPair[1] = gsl_quat_float_alloc();

  gsl_vector_float_memcpy(rot->value.quat_invPair[0], pQ);
  gsl_vector_float_memcpy(rot->value.quat_invPair[1], pQ);

  

  return gsl_quat_float_inv(rot->value.quat_invPair[1]);
}

/**
 *
 * Private function definitions
 *
 */