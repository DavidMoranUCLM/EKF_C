#include "jacobian.h"




void jacobian(gsl_matrix_float *pJacMatrix, diff_function pFunc, gsl_vector_float *pVec0, float delta){
    gsl_vector_float* pVecDelta = gsl_vector_float_alloc(pVec0->size);
    gsl_vector_float* pVecForwardResult = gsl_vector_float_alloc(pVec0->size);
    gsl_vector_float* pVecBackwardResult = gsl_vector_float_alloc(pVec0->size);

    for(size_t i=0;i<pVec0->size;i++){
        gsl_vector_float_memcpy(pVecDelta, pVec0);
        gsl_vector_float_set(pVecDelta, i, gsl_vector_float_get(pVec0,i)+delta);

        pFunc(pVecDelta, pVecForwardResult);

        gsl_vector_float_memcpy(pVecDelta, pVec0);
        gsl_vector_float_set(pVecDelta, i, gsl_vector_float_get(pVec0,i)-delta);

        pFunc(pVecDelta, pVecBackwardResult);

        gsl_vector_float_sub(pVecForwardResult,pVecBackwardResult);
        gsl_matrix_float_set_col(pJacMatrix, i, pVecForwardResult);
    }

    gsl_matrix_float_scale(pJacMatrix, 0.5/delta);

    gsl_vector_float_free(pVecDelta);
    gsl_vector_float_free(pVecForwardResult);
    gsl_vector_float_free(pVecBackwardResult);
}