#include "unity.h"

#include <math.h>
#include <stdlib.h>

#include "gsl/gsl_math.h"
#include "gsl_quaternion_float.h"

#define PROD_NORM_ITER 10000
#define FLOAT_ERROR 1e-6

#define TEST_ASSERT_QUAT(q1, q2) \
  TEST_ASSERT_EQUAL_FLOAT_ARRAY((q1->data), (q2->data), (4))

void setUp(void) {}
void tearDown(void) {}

void suiteSetUp(void) {}
//int suiteTearDown(int num_failures) {}

void resetTest(void) {}
void verifyTest(void) {}

void randQuat(gsl_quat_float* q) {
  for (size_t i = 0; i < q->size; i++) {
    q->data[i] = ((float)rand() * 100) / (RAND_MAX);
  }
}

void isNormalized(gsl_quat_float* q) {
  float norm = gsl_quat_float_norm(q);
  TEST_ASSERT_EQUAL_FLOAT(1, norm);
}

void testQuatAlloc(void) {
  gsl_quat_float* q = gsl_quat_float_alloc();

  TEST_ASSERT_NOT_NULL(q);

  gsl_quat_float_free(q);
}

void testQuatCalloc(void) {
  gsl_quat_float* q = gsl_quat_float_calloc();

  TEST_ASSERT_NOT_NULL(q);
  for (size_t i = 0; i < q->size; i++) {
    TEST_ASSERT(q->data[i] == 0);
  }

  gsl_quat_float_free(q);
}

void testQuatNorm(void) {
  gsl_quat_float* q = gsl_quat_float_calloc();
  for (int i = 0; i < 100; i++) {
    float normRef = 0;
    float norm;

    randQuat(q);

    for (size_t i = 0; i < q->size; i++) {
      normRef += q->data[i] * q->data[i];
    }

    normRef = sqrtf(normRef);

    norm = gsl_quat_float_norm(q);

    TEST_ASSERT_FLOAT_ARRAY_WITHIN(FLOAT_ERROR, &normRef, &norm, 1);
  }

  gsl_quat_float_free(q);
}

void testQuatNormalize(void) {
  gsl_quat_float* q = gsl_quat_float_alloc();

  for (size_t i = 0; i < q->size; i++) {
    q->data[i] = ((float)rand() * 100) / (RAND_MAX);
  }

  gsl_quat_float_normalize(q);

  isNormalized(q);
}

void testQuatConjugate(void) {
  gsl_quat_float* pQ1 = gsl_quat_float_alloc();
  gsl_quat_float* pQ1ConjRef = gsl_quat_float_alloc();

  randQuat(pQ1);

  pQ1ConjRef->data[0] = pQ1->data[0];
  pQ1ConjRef->data[1] = -pQ1->data[1];
  pQ1ConjRef->data[2] = -pQ1->data[2];
  pQ1ConjRef->data[3] = -pQ1->data[3];

  gsl_quat_float_conjugate(pQ1);

  TEST_ASSERT_QUAT(pQ1ConjRef, pQ1);

  gsl_quat_float_free(pQ1);
  gsl_quat_float_free(pQ1ConjRef);
}

void testQuatCopy(void) {
  gsl_quat_float* q1 = gsl_quat_float_alloc();
  gsl_quat_float* q2 = gsl_quat_float_alloc();
  gsl_quat_float* q3 = gsl_quat_float_alloc();

  randQuat(q1);
  randQuat(q2);
  randQuat(q3);

  gsl_quat_float_copy(q1, q2);
  gsl_quat_float_copy(q1, q3);

  TEST_ASSERT_QUAT(q2, q3);
}

void testQuatProd(void) {
  gsl_quat_float* q1 = gsl_quat_float_alloc();
  gsl_quat_float* q2 = gsl_quat_float_alloc();
  gsl_quat_float* q3Ref = gsl_quat_float_alloc();

  for (uint8_t i = 0; i < q1->size; i++) {
    q1->data[i] = i + 1;
    q2->data[i] = i + 1;
  }

  q3Ref->data[0] = -28;
  q3Ref->data[1] = 4;
  q3Ref->data[2] = 6;
  q3Ref->data[3] = 8;

  gsl_quat_float_product(q1, q2);

  TEST_ASSERT_QUAT(q3Ref, q1);

  gsl_quat_float_free(q1);
  gsl_quat_float_free(q2);
  gsl_quat_float_free(q3Ref);
}

void testQuatProdNorm(void) {
  gsl_quat_float* q1 = gsl_quat_float_alloc();
  gsl_quat_float* q2 = gsl_quat_float_alloc();

  randQuat(q1);
  randQuat(q2);

  gsl_quat_float_normalize(q1);
  gsl_quat_float_normalize(q2);

  for (int i = 0; i < PROD_NORM_ITER; i++) {
    gsl_quat_float_product(q1, q2);
    gsl_quat_float_normalize(q1);
    isNormalized(q1);
    randQuat(q1);
    randQuat(q2);
  }

  gsl_quat_float_free(q1);
  gsl_quat_float_free(q2);
}

void testQuatFromAxis(void) {
  gsl_quat_float* pQ1 = gsl_quat_float_calloc();
  gsl_quat_float* pQ1Ref = gsl_quat_float_calloc();

  gsl_vector_float* pAxis = gsl_vector_float_calloc(3);

  float angle = M_PI / 4;
  pAxis->data[0] = 0;
  pAxis->data[1] = 0;
  pAxis->data[2] = 1;

  pQ1Ref->data[0] = 0.9238795;
  pQ1Ref->data[1] = 0;
  pQ1Ref->data[2] = 0;
  pQ1Ref->data[3] = 0.3826835;

  TEST_ASSERT(gsl_quat_float_fromAxis(pAxis, angle, pQ1)==0);

  TEST_ASSERT_QUAT(pQ1Ref, pQ1);

  pAxis->data[0] = 0;
  pAxis->data[1] = 0;
  pAxis->data[2] = 0;

  TEST_ASSERT_FALSE(gsl_quat_float_fromAxis(pAxis, angle, pQ1)==0);

  gsl_quat_float_free(pQ1);
  gsl_quat_float_free(pQ1Ref);
  gsl_quat_float_free(pAxis);
}

void testQuatFromVector(void) {
  gsl_quat_float* pQ1 = gsl_quat_float_alloc();
  gsl_vector_float* pV1 = gsl_vector_float_alloc(3);

  randQuat(pV1);

  TEST_ASSERT(gsl_quat_float_fromVector(pV1,pQ1)==0);

  TEST_ASSERT_EQUAL_FLOAT(.0f, pQ1->data[0]);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY(pV1->data, pQ1->data + 1, 3);

  free(pV1);

  pV1 = gsl_vector_float_alloc(4);
  randQuat(pV1);

  TEST_ASSERT_FALSE(gsl_quat_float_fromVector(pV1,pQ1)==0);
}

void testQuatToMatrix(void){
  gsl_quat_float *pQ = gsl_quat_float_alloc();
  randQuat(pQ);

  gsl_matrix_float *pMatrix = gsl_matrix_float_alloc(4,4);
  gsl_quat_float_toMatrix(pQ, pMatrix);

  TEST_ASSERT_EQUAL_FLOAT(gsl_matrix_float_get(pMatrix,0,0), gsl_quat_float_get(pQ,0));
  TEST_ASSERT_EQUAL_FLOAT(gsl_matrix_float_get(pMatrix,0,1), -gsl_quat_float_get(pQ,1));
  TEST_ASSERT_EQUAL_FLOAT(gsl_matrix_float_get(pMatrix,0,2), -gsl_quat_float_get(pQ,2));
  TEST_ASSERT_EQUAL_FLOAT(gsl_matrix_float_get(pMatrix,0,3), -gsl_quat_float_get(pQ,3));
  TEST_ASSERT_EQUAL_FLOAT(gsl_matrix_float_get(pMatrix,1,1), gsl_quat_float_get(pQ,0));
  TEST_ASSERT_EQUAL_FLOAT(gsl_matrix_float_get(pMatrix,1,2), gsl_quat_float_get(pQ,3));
  TEST_ASSERT_EQUAL_FLOAT(gsl_matrix_float_get(pMatrix,1,3), -gsl_quat_float_get(pQ,2));
  TEST_ASSERT_EQUAL_FLOAT(gsl_matrix_float_get(pMatrix,2,2), gsl_quat_float_get(pQ,0));
  TEST_ASSERT_EQUAL_FLOAT(gsl_matrix_float_get(pMatrix,2,3), gsl_quat_float_get(pQ,1));
  TEST_ASSERT_EQUAL_FLOAT(gsl_matrix_float_get(pMatrix,3,3), gsl_quat_float_get(pQ,0));

  for(uint8_t i=0; i<4; i++){
    for(uint8_t j=0; j<4; j++){
      TEST_ASSERT_EQUAL_FLOAT(gsl_matrix_float_get(pMatrix,i,j), gsl_matrix_float_get(pMatrix,i,j));
    }
  }
}

void testQuatToRotMatrix(void) {
  // Create a quaternion for 90° rotation about Z axis: [cos(pi/4), 0, 0, sin(pi/4)]
  gsl_quat_float* q = gsl_quat_float_alloc();
  q->data[0] = 0.70710678f; // cos(pi/4)
  q->data[1] = 0.0f;
  q->data[2] = 0.0f;
  q->data[3] = 0.70710678f; // sin(pi/4)

  // Allocate a 3x3 rotation matrix
  gsl_matrix_float* rotMat = gsl_matrix_float_alloc(3, 3);
  TEST_ASSERT(gsl_quat_float_toRotMatrix(q, rotMat) == 0);

  // Expected rotation matrix for 90° about Z:
  // [ 0, -1, 0 ]
  // [ 1,  0, 0 ]
  // [ 0,  0, 1 ]
  float expected[9] = {
      0.f, -1.f, 0.f,
      1.f,  0.f, 0.f,
      0.f,  0.f, 1.f
  };

  for (size_t i = 0; i < 3; i++) {
    for (size_t j = 0; j < 3; j++) {
      float val = gsl_matrix_float_get(rotMat, i, j);
      TEST_ASSERT_FLOAT_WITHIN(FLOAT_ERROR, expected[i * 3 + j], val);
    }
  }

  gsl_matrix_float_free(rotMat);
  gsl_quat_float_free(q);
}

void testQuatToRotMatrixIdentity(void) {
  // Expected identity matrix 3x3
  float identity[9] = {
    1.f, 0.f, 0.f,
    0.f, 1.f, 0.f,
    0.f, 0.f, 1.f
  };

  for (int sign = 1; sign >= -1; sign -= 2) {
    gsl_quat_float* q = gsl_quat_float_alloc();
    q->data[0] = (float)sign;  // +1 or -1
    q->data[1] = 0.f;
    q->data[2] = 0.f;
    q->data[3] = 0.f;

    gsl_matrix_float* rotMat = gsl_matrix_float_alloc(3, 3);
    TEST_ASSERT(gsl_quat_float_toRotMatrix(q, rotMat) == 0);

    for (size_t i = 0; i < 3; i++) {
      for (size_t j = 0; j < 3; j++) {
        float val = gsl_matrix_float_get(rotMat, i, j);
        TEST_ASSERT_FLOAT_WITHIN(FLOAT_ERROR, identity[i * 3 + j], val);
      }
    }
    gsl_matrix_float_free(rotMat);
    gsl_quat_float_free(q);
  }
}

int main(void) {
  UNITY_BEGIN();
  RUN_TEST(testQuatAlloc);
  RUN_TEST(testQuatCalloc);
  RUN_TEST(testQuatNorm);
  RUN_TEST(testQuatNormalize);
  RUN_TEST(testQuatConjugate);
  RUN_TEST(testQuatCopy);
  RUN_TEST(testQuatProd);
  RUN_TEST(testQuatProdNorm);
  RUN_TEST(testQuatFromAxis);
  RUN_TEST(testQuatFromVector);
  RUN_TEST(testQuatToMatrix);
  RUN_TEST(testQuatToRotMatrix);
  RUN_TEST(testQuatToRotMatrixIdentity);
  return UNITY_END();
}