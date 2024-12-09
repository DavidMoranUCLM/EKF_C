#include "math_utils.h"


int skewSymFromVector(gsl_vector_float *pV, gsl_matrix_float *pM){
  if ((pV->size != pM->size1) || (pV->size != pM->size2)){
    return -1;
  }

  if (pV->size != 3){
    return -1;
  }

  gsl_matrix_float_set_zero(pM);
    
  gsl_matrix_float_set(pM,0,1, -gsl_vector_float_get(pV,2));
  gsl_matrix_float_set(pM,0,2, gsl_vector_float_get(pV,1));
  gsl_matrix_float_set(pM,1,2, -gsl_vector_float_get(pV,0));

  gsl_matrix_float *pMT = gsl_matrix_float_alloc(3,3);
  gsl_matrix_float_transpose_memcpy(pMT,pM);
  gsl_matrix_float_scale(pMT,-1.f);

  gsl_matrix_float_add(pM,pMT);
  gsl_matrix_float_free(pMT);

  return 0;

}

void cross_product(const gsl_vector_float *a, const gsl_vector_float *b,
                   gsl_vector_float *c) {
  if (a->size != 3 || b->size != 3 || c->size != 3) {
    fprintf(stderr, "Error: Vectors must have size 3.\n");
    return;
  }

  float c1 = gsl_vector_float_get(a, 1) * gsl_vector_float_get(b, 2) -
             gsl_vector_float_get(a, 2) * gsl_vector_float_get(b, 1);
  float c2 = gsl_vector_float_get(a, 2) * gsl_vector_float_get(b, 0) -
             gsl_vector_float_get(a, 0) * gsl_vector_float_get(b, 2);
  float c3 = gsl_vector_float_get(a, 0) * gsl_vector_float_get(b, 1) -
             gsl_vector_float_get(a, 1) * gsl_vector_float_get(b, 0);

  gsl_vector_float_set(c, 0, c1);
  gsl_vector_float_set(c, 1, c2);
  gsl_vector_float_set(c, 2, c3);
}