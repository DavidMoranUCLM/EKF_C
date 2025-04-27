
#include "mag_correction.h"

#include "gsl/gsl_blas.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_math.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_matrix_float.h"
#include "gsl/gsl_vector.h"
#include "gsl/gsl_vector_float.h"
#include "gsl_nan.h"
#include "math.h"
#include "math_utils.h"
#include "utils_rotmatYawCorrectionJac.h"

void qCorrectMag(gsl_vector_float *state, const gsl_vector_float *mag) {
  // Extract quaternion components
  float w = gsl_vector_float_get(state, 0);
  float x = gsl_vector_float_get(state, 1);
  float y = gsl_vector_float_get(state, 2);
  float z = gsl_vector_float_get(state, 3);

  // Construct rotation matrix R from quaternion
  gsl_matrix_float *R = gsl_matrix_float_alloc(3, 3);
  gsl_matrix_float_set(R, 0, 0, 1 - 2 * (y * y + z * z));
  gsl_matrix_float_set(R, 0, 1, 2 * (x * y - z * w));
  gsl_matrix_float_set(R, 0, 2, 2 * (x * z + y * w));
  
  gsl_matrix_float_set(R, 1, 0, 2 * (x * y + z * w));
  gsl_matrix_float_set(R, 1, 1, 1 - 2 * (x * x + z * z));
  gsl_matrix_float_set(R, 1, 2, 2 * (y * z - x * w));
  
  gsl_matrix_float_set(R, 2, 0, 2 * (x * z - y * w));
  gsl_matrix_float_set(R, 2, 1, 2 * (y * z + x * w));
  gsl_matrix_float_set(R, 2, 2, 1 - 2 * (x * x + y * y));

  // Compute yaw from rotation matrix
  float R11 = gsl_matrix_float_get(R, 0, 0);
  float R21 = gsl_matrix_float_get(R, 1, 0);
  float yaw = atan2f(R21, R11);

  // Create unyaw rotation matrix (around Z-axis by -yaw)
  gsl_matrix_float *R_unyaw = gsl_matrix_float_alloc(3, 3);
  float cy = cosf(-yaw), sy = sinf(-yaw);
  gsl_matrix_float_set(R_unyaw, 0, 0, cy);
  gsl_matrix_float_set(R_unyaw, 0, 1, -sy);
  gsl_matrix_float_set(R_unyaw, 0, 2, 0);
  gsl_matrix_float_set(R_unyaw, 1, 0, sy);
  gsl_matrix_float_set(R_unyaw, 1, 1, cy);
  gsl_matrix_float_set(R_unyaw, 1, 2, 0);
  gsl_matrix_float_set(R_unyaw, 2, 0, 0);
  gsl_matrix_float_set(R_unyaw, 2, 1, 0);
  gsl_matrix_float_set(R_unyaw, 2, 2, 1);

  // Compute R_tilt = R_unyaw * R
  gsl_matrix_float *R_tilt = gsl_matrix_float_alloc(3, 3);
  gsl_blas_sgemm(CblasNoTrans, CblasNoTrans, 1.0, R_unyaw, R, 0.0, R_tilt);

  gsl_vector_float *mag_corrected = gsl_vector_float_alloc(3);
  gsl_blas_sgemv(CblasNoTrans, 1.0, R_tilt, mag, 0.0, mag_corrected);

  // Compute new yaw from corrected magnetometer values
  float mc_x = gsl_vector_float_get(mag_corrected, 0);
  float mc_y = gsl_vector_float_get(mag_corrected, 1);
  float new_yaw = -atan2f(mc_y, mc_x) + M_PI_2;

  // Create yaw correction matrix (around Z-axis by new_yaw)
  gsl_matrix_float *R_yaw = gsl_matrix_float_alloc(3, 3);
  float cny = cosf(new_yaw), sny = sinf(new_yaw);
  gsl_matrix_float_set(R_yaw, 0, 0, cny);
  gsl_matrix_float_set(R_yaw, 0, 1, -sny);
  gsl_matrix_float_set(R_yaw, 0, 2, 0);
  gsl_matrix_float_set(R_yaw, 1, 0, sny);
  gsl_matrix_float_set(R_yaw, 1, 1, cny);
  gsl_matrix_float_set(R_yaw, 1, 2, 0);
  gsl_matrix_float_set(R_yaw, 2, 0, 0);
  gsl_matrix_float_set(R_yaw, 2, 1, 0);
  gsl_matrix_float_set(R_yaw, 2, 2, 1);

  // Compute corrected rotation matrix R_corrected = R_yaw * R_tilt
  gsl_matrix_float *R_corrected = gsl_matrix_float_alloc(3, 3);
  gsl_blas_sgemm(CblasNoTrans, CblasNoTrans, 1.0, R_yaw, R_tilt, 0.0,
                 R_corrected);

  // Convert R_corrected to quaternion
  float qw, qx, qy, qz;
  float tr = gsl_matrix_float_get(R_corrected, 0, 0) +
             gsl_matrix_float_get(R_corrected, 1, 1) +
             gsl_matrix_float_get(R_corrected, 2, 2);
  float R32 = gsl_matrix_float_get(R_corrected, 2, 1),
        R23 = gsl_matrix_float_get(R_corrected, 1, 2);
  float R13 = gsl_matrix_float_get(R_corrected, 0, 2),
        R31 = gsl_matrix_float_get(R_corrected, 2, 0);
  float R12 = gsl_matrix_float_get(R_corrected, 0, 1);
  R21 = gsl_matrix_float_get(R_corrected, 1, 0);

  if (tr > 0) {
    qw = 0.5f * sqrtf(1 + tr);
    float denom = 4 * qw;
    qx = (R32 - R23) / denom;
    qy = (R13 - R31) / denom;
    qz = (R21 - R12) / denom;
  } else {
    float d0 = gsl_matrix_float_get(R_corrected, 0, 0);
    float d1 = gsl_matrix_float_get(R_corrected, 1, 1);
    float d2 = gsl_matrix_float_get(R_corrected, 2, 2);

    if (d0 >= d1 && d0 >= d2) {
      float s = 2 * sqrtf(1 + d0 - d1 - d2);
      qw = (R32 - R23) / s;
      qx = 0.25 * s;
      qy = (R12 + R21) / s;
      qz = (R13 + R31) / s;
    } else if (d1 >= d2) {
      float s = 2 * sqrtf(1 + d1 - d0 - d2);
      qw = (R13 - R31) / s;
      qx = (R12 + R21) / s;
      qy = 0.25 * s;
      qz = (R23 + R32) / s;
    } else {
      float s = 2 * sqrtf(1 + d2 - d0 - d1);
      qw = (R21 - R12) / s;
      qx = (R13 + R31) / s;
      qy = (R23 + R32) / s;
      qz = 0.25 * s;
    }
  }

  gsl_vector_float_set(state, 0, qw);
  gsl_vector_float_set(state, 1, qx);
  gsl_vector_float_set(state, 2, qy);
  gsl_vector_float_set(state, 3, qz);

  // Free allocated GSL objects
  gsl_matrix_float_free(R);
  gsl_matrix_float_free(R_unyaw);
  gsl_matrix_float_free(R_tilt);
  gsl_vector_float_free(mag_corrected);
  gsl_matrix_float_free(R_yaw);
  gsl_matrix_float_free(R_corrected);
}

void PCorrectMag(const gsl_vector_float *x, const gsl_vector_float *mag,
                 const gsl_vector_float *sigma_mag, gsl_matrix_float *P) {
  gsl_matrix_float *J = gsl_matrix_float_alloc(7, 3);
  gsl_matrix_float *tmp4_3 = gsl_matrix_float_alloc(4, 3);
  gsl_matrix_float_view J_4_3 = gsl_matrix_float_submatrix(J, 0, 0, 4, 3);
  gsl_matrix_float *sigma = gsl_matrix_float_calloc(3, 3);

  // float J1[3][7];
  utils_rotmatYawCorrectionJac(x->data, mag->data, J->data);
  for (int i = 0; i < 7; ++i) {
    for (int j = 0; j < 3; ++j) {
      if (gsl_isnan(gsl_matrix_float_get(J, i, j))) {
        gsl_matrix_float_set(J, i, j, 10000.f);
      }
    }
  }

  gsl_matrix_float_set(sigma, 0, 0, sigma_mag->data[0]);
  gsl_matrix_float_set(sigma, 1, 1, sigma_mag->data[1]);
  gsl_matrix_float_set(sigma, 2, 2, sigma_mag->data[2]);

  gsl_blas_sgemm(CblasNoTrans, CblasNoTrans, 1.0, &J_4_3.matrix, sigma, 0.0,
                 tmp4_3);
  gsl_blas_sgemm(CblasNoTrans, CblasTrans, 1.0, tmp4_3, &J_4_3.matrix, 0.0, P);

  gsl_matrix_float_free(J);
  gsl_matrix_float_free(tmp4_3);
  gsl_matrix_float_free(sigma);
  return;
}

void correctMag(gsl_matrix_float *P1, gsl_vector_float *x1,
                const gsl_vector_float *mag,
                const gsl_vector_float *mag_sigma) {
  gsl_matrix_float *P2 = gsl_matrix_float_alloc(4, 4);
  gsl_matrix_float *P3 = gsl_matrix_float_alloc(4, 4);

  gsl_vector_float *x2 = gsl_vector_float_alloc(4);
  gsl_vector_float *x3 = gsl_vector_float_alloc(4);

  gsl_vector_float_memcpy(x2, x1);
  qCorrectMag(x2, mag);
  PCorrectMag(x1, mag, mag_sigma, P2);

  normal_dist_intersection(x1, x2, x3, P1, P2, P3);

  gsl_vector_float_memcpy(x1, x3);
  gsl_matrix_float_memcpy(P1, P3);

  gsl_matrix_float_free(P2);
  gsl_matrix_float_free(P3);
  gsl_vector_float_free(x2);
  gsl_vector_float_free(x3);
  return;
}