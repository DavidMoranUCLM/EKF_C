#include <stdio.h>
#include "mag_correction.h"
#include "gsl/gsl_matrix_float.h"
#include "gsl/gsl_vector_float.h"
#include "unity.h"
#include <time.h>
#include <stdlib.h>

#define MAG_SIGMA 0.01f*0.01f

gsl_vector_float *sigma_mag = NULL;

// Helper to allocate and initialize a quaternion
static gsl_vector_float* create_quaternion(float w, float x, float y, float z) {
  gsl_vector_float *q = gsl_vector_float_alloc(4);
  gsl_vector_float_set(q, 0, w);
  gsl_vector_float_set(q, 1, x);
  gsl_vector_float_set(q, 2, y);
  gsl_vector_float_set(q, 3, z);
  return q;
}

// Helper to allocate and initialize a 3x3 identity matrix
static gsl_matrix_float* create_identity_matrix() {
  gsl_matrix_float *P = gsl_matrix_float_alloc(4, 4);
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      gsl_matrix_float_set(P, i, j, (i == j) ? 1.0f : 0.0f);
  return P;
}

// Helper to allocate and initialize a vector
static gsl_vector_float* create_vector(float x, float y, float z) {
  gsl_vector_float *v = gsl_vector_float_alloc(3);
  gsl_vector_float_set(v, 0, x);
  gsl_vector_float_set(v, 1, y);
  gsl_vector_float_set(v, 2, z);
  return v;
}

void setUp(void) {}
void tearDown(void) {}

void test_mag_correct_identity_quaternion(void) {
  gsl_vector_float *q = create_quaternion(1.0f, 0.0f, 0.0f, 0.0f);
  gsl_matrix_float *P = create_identity_matrix();
  gsl_vector_float *mag = create_vector(0.5f, 0.5f, 0.0f);

  correctMag(P, q, mag, sigma_mag);


  gsl_vector_float_free(q);
  gsl_matrix_float_free(P);
  gsl_vector_float_free(mag);
}

void test_mag_correct_non_identity_quaternion(void) {
  gsl_vector_float *q = create_quaternion(0.7071f, 0.7071f, 0.0f, 0.0f);
  gsl_matrix_float *P = create_identity_matrix();
  gsl_vector_float *mag = create_vector(0.5f, 0.5f, 0.0f);

  correctMag(P, q, mag, sigma_mag);

  float norm = 0.0f;
  for (int i = 0; i < 4; ++i)
    norm += gsl_vector_float_get(q, i) * gsl_vector_float_get(q, i);
  TEST_ASSERT_FLOAT_WITHIN(1e-3, 1.0f, sqrtf(norm));

  gsl_vector_float_free(q);
  gsl_matrix_float_free(P);
  gsl_vector_float_free(mag);
}

void test_mag_correct_zero_magnetometer(void) {
  gsl_vector_float *q = create_quaternion(1.0f, 0.0f, 0.0f, 0.0f);
  gsl_matrix_float *P = create_identity_matrix();
  gsl_vector_float *mag = create_vector(0.0f, 0.0f, 0.0f);

  correctMag(P, q, mag, sigma_mag);

  // Quaternion should remain normalized
  float norm = 0.0f;
  for (int i = 0; i < 4; ++i)
    norm += gsl_vector_float_get(q, i) * gsl_vector_float_get(q, i);
  TEST_ASSERT_FLOAT_WITHIN(1e-3, 1.0f, sqrtf(norm));
  gsl_vector_float_free(q);
  gsl_matrix_float_free(P);
  gsl_vector_float_free(mag);
}

void test_mag_correct_covariance_update(void) {
  gsl_vector_float *q = create_quaternion(1.0f, 0.0f, 0.0f, 0.0f);
  gsl_matrix_float *P = create_identity_matrix();
  gsl_vector_float *mag = create_vector(1.0f, 0.0f, 0.0f);

  correctMag(P, q, mag, sigma_mag);

  // Covariance diagonal should remain positive
  for (int i = 0; i < 3; ++i)
    TEST_ASSERT_TRUE(gsl_matrix_float_get(P, i, i) > 0.0f);

  gsl_vector_float_free(q);
  gsl_matrix_float_free(P);
  gsl_vector_float_free(mag);
}

void test_mag_correct_timing(void) {
  const int N = 10000;
  gsl_vector_float *q = gsl_vector_float_alloc(4);
  gsl_matrix_float *P = gsl_matrix_float_alloc(4, 4);
  gsl_vector_float *mag = gsl_vector_float_alloc(3);
  gsl_matrix_float *A = gsl_matrix_float_alloc(4, 4);
  srand((unsigned)time(NULL));

  clock_t start = clock();
  for (int i = 0; i < N; ++i) {
    // Random quaternion
    float qnorm = 0.0f;
    for (int j = 0; j < 4; ++j) {
      float val = ((float)rand() / RAND_MAX) * 2.0f - 1.0f;
      gsl_vector_float_set(q, j, val);
      qnorm += val * val;
    }
    qnorm = sqrtf(qnorm);
    for (int j = 0; j < 4; ++j)
      gsl_vector_float_set(q, j, gsl_vector_float_get(q, j) / qnorm);
    // Random positive definite 4x4 matrix P = A * A^T
    for (int r = 0; r < 4; ++r)
      for (int c = 0; c < 4; ++c)
        gsl_matrix_float_set(A, r, c, ((float)rand() / RAND_MAX) * 2.0f - 1.0f);
    for (int r = 0; r < 4; ++r) {
      for (int c = 0; c < 4; ++c) {
        float sum = 0.0f;
        for (int k = 0; k < 4; ++k)
          sum += gsl_matrix_float_get(A, r, k) * gsl_matrix_float_get(A, c, k);
        gsl_matrix_float_set(P, r, c, sum);
      }
    }
    // Random magnetometer
    float mnorm = 0.0f;
    for (int j = 0; j < 3; ++j) {
      float val = ((float)rand() / RAND_MAX) * 2.0f - 1.0f;
      gsl_vector_float_set(mag, j, val);
      mnorm += val * val;
    }
    mnorm = sqrtf(mnorm);
    if (mnorm > 1e-6f) {
      for (int j = 0; j < 3; ++j)
        gsl_vector_float_set(mag, j, gsl_vector_float_get(mag, j) / mnorm);
    }
    correctMag(P, q, mag, sigma_mag);
  }
  clock_t end = clock();
  double elapsed = (double)(end - start) / CLOCKS_PER_SEC;
  printf("Time for %d executions of correctMag: %f seconds\n", N, elapsed);

  gsl_vector_float_free(q);
  gsl_matrix_float_free(P);
  gsl_vector_float_free(mag);
  gsl_matrix_float_free(A);
}

void test_normal_fusion(void) {
  // Test for normal_fusion (normal_dist_intersection)
  gsl_vector_float *v1 = gsl_vector_float_alloc(4);
  gsl_vector_float *v2 = gsl_vector_float_alloc(4);
  gsl_vector_float *v3 = gsl_vector_float_alloc(4);
  gsl_matrix_float *P1 = gsl_matrix_float_alloc(4, 4);
  gsl_matrix_float *P2 = gsl_matrix_float_alloc(4, 4);
  gsl_matrix_float *P3 = gsl_matrix_float_alloc(4, 4);

  // Set v1 and v2 to known values
  for (int i = 0; i < 4; ++i) {
    gsl_vector_float_set(v1, i, 1.0f + i);
    gsl_vector_float_set(v2, i, 2.0f + i);
    for (int j = 0; j < 4; ++j) {
      gsl_matrix_float_set(P1, i, j, (i == j) ? 1.0f +j: 0.0f);
      gsl_matrix_float_set(P2, i, j, (i == j) ? 2.0f +j: 0.0f);
    }
  }

  // Call the function under test
  int status = normal_dist_intersection(v1, v2, v3, P1, P2, P3);
  TEST_ASSERT_EQUAL_INT(0, status);

  // Check that v3 and P3 are finite and P3 diagonal is positive
  for (int i = 0; i < 4; ++i) {
    TEST_ASSERT_TRUE(isfinite(gsl_vector_float_get(v3, i)));
    TEST_ASSERT_TRUE(gsl_matrix_float_get(P3, i, i) > 0.0f);
  }

  gsl_vector_float_free(v1);
  gsl_vector_float_free(v2);
  gsl_vector_float_free(v3);
  gsl_matrix_float_free(P1);
  gsl_matrix_float_free(P2);
  gsl_matrix_float_free(P3);
}

int main(void) {
  sigma_mag = (gsl_vector_float*)gsl_vector_float_alloc(3);
  gsl_vector_float_set_all(sigma_mag, MAG_SIGMA);

  UNITY_BEGIN();
  RUN_TEST(test_mag_correct_identity_quaternion);
  RUN_TEST(test_mag_correct_non_identity_quaternion);
  RUN_TEST(test_mag_correct_zero_magnetometer);
  RUN_TEST(test_mag_correct_covariance_update);
  RUN_TEST(test_mag_correct_timing);
  RUN_TEST(test_normal_fusion);
  
  return UNITY_END();
}