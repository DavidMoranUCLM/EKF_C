#include "jacobian.h"
#include "math.h"
#include "unity.h"

void setUp(void) {}
void tearDown(void) {}

void suiteSetUp(void) {}
//int suiteTearDown(int num_failures) {}

void resetTest(void) {}
void verifyTest(void) {}

//f(x1,x2,x3) = [x1^2, x2^2, x3^2]
void func(gsl_vector_float *in, gsl_vector_float *out);

void test_jacobian1(){
    float elem1 = 2;
    float elem2 = 5;
    float elem3 = 2;
    gsl_vector_float *pCenterVec = gsl_vector_float_alloc(3);
    gsl_vector_float_set(pCenterVec, 0, elem1);
    gsl_vector_float_set(pCenterVec, 1, elem2);
    gsl_vector_float_set(pCenterVec, 2, elem3);

    gsl_matrix_float *pJac = gsl_matrix_float_alloc(3,3);

    jacobian(pJac,func,pCenterVec, 5e-3F);

    float expectedResult[3][3] = {{elem1*2, 0, 0},
                                  {0, elem2*2, 0},
                                  {0, 0, elem3*2}};
    
    TEST_ASSERT_FLOAT_ARRAY_WITHIN(1e-2F,expectedResult, pJac->data, 9);
}

int main(){
    UNITY_BEGIN();
    RUN_TEST(test_jacobian1);
    return UNITY_END();
}


void func(gsl_vector_float *in, gsl_vector_float *out){
    float elem1, elem2, elem3;

    elem1 = pow(gsl_vector_float_get(in,0),2);
    elem2 = pow(gsl_vector_float_get(in,1),2);
    elem3 = pow(gsl_vector_float_get(in,2),2);

    gsl_vector_float_set(out, 0, elem1);
    gsl_vector_float_set(out, 1, elem2);
    gsl_vector_float_set(out, 2, elem3);

}