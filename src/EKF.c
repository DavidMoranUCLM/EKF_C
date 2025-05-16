#include "EKF.h"

#include "EKF_const.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_const.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_math.h"
#include "gsl/gsl_matrix_float.h"
#include "gsl/gsl_vector_float.h"
#include "gsl_quaternion_float.h"
#include "mag_correction.h"
#include "math.h"
#include "math_utils.h"
#include "stdio.h"

/**
 *  Private Definitions
 */

#define PRINT_MATRIX(m) gsl_matrix_float_fprintf(stdout, m, "%f")

/**
 * Private Functions Declarations
 */

void initWorkSpace(EKF_ctx_t* ctx);
void deInitWorkSpace(EKF_ctx_t* ctx);
void ekfUpdate(EKF_ctx_t* ctx, const measures_t* measures,
               const float currentTime);

void ekfEstimate(EKF_ctx_t* ctx);
void qEst(EKF_ctx_t* ctx);
void PEst(EKF_ctx_t* ctx);
void getF(EKF_ctx_t* ctx);
void getQ(EKF_ctx_t* ctx);
void getW(EKF_ctx_t* ctx);

void ekfCorrect(EKF_ctx_t* ctx);
void PCorrect(EKF_ctx_t* ctx);
void get_h(EKF_ctx_t* ctx);
int getH(EKF_ctx_t* ctx);
int getR(EKF_ctx_t* ctx);
int getS(EKF_ctx_t* ctx);
int getK(EKF_ctx_t* ctx);

void qEstPrimitive(const gsl_vector_float* velAng, float deltaT,
                   const gsl_quat_float* qPrev, gsl_quat_float* qEst,
                   gsl_quat_float* tmpQuat, gsl_matrix_float* qVelAngMat,
                   gsl_matrix_float* q1Mat);
void PEstPrimitive(const gsl_matrix_float* PPrev, const gsl_matrix_float* F,
                   const gsl_matrix_float* Q, gsl_matrix_float* PEst,
                   gsl_matrix_float* tmp4x4);
void PCorrectPrimitive(const gsl_matrix_float* P, const gsl_matrix_float* K,
                       const gsl_matrix_float* H, const gsl_matrix_float* R,
                       gsl_matrix_float* PCorrect, gsl_matrix_float* tmp4x4,
                       gsl_matrix_float* tmp3x4, gsl_matrix_float* I4);
void qCorrectPrimitive(const gsl_quat_float* q, const gsl_matrix_float* K,
                       const gsl_vector_float* h, gsl_quat_float* qCorrect,
                       gsl_vector_float* tmp);
void getKPrimitive(const gsl_matrix_float* P, const gsl_matrix_float* H,
                   const gsl_matrix_float* invS, gsl_matrix_float* K,
                   gsl_matrix_float* tmp4x3);

void getHPrimitive(const gsl_quat_float* q, const gsl_vector_float* acc,
                   gsl_matrix_float* H, gsl_vector_float* pQv,
                   gsl_matrix_float* pM2, gsl_vector_float* pV1,
                   gsl_vector_float* pV2);
void getRPrimitive(gsl_matrix_float* R);
void getSPrimitive(const gsl_matrix_float* H, const gsl_matrix_float* P,
                   const gsl_matrix_float* R, gsl_matrix_float* S,
                   gsl_matrix_float* tmp4x3);
void getFPrimitive(const gsl_vector_float* velAng, float deltaT,
                   gsl_vector_float* prevState, gsl_quat_float* tmpQuat,
                   gsl_matrix_float* F, const gsl_matrix_float* I4);
void getWPrimitive(const gsl_quat_float* qPrev, float deltaT,
                   gsl_matrix_float* W);
void get_hPrimitive(const gsl_vector_float* state_est,
                    const gsl_vector_float* horizonRefG, gsl_vector_float* h);

void invertMatrixFloat(EKF_ctx_t* ctx, const gsl_matrix_float* S,
                       gsl_matrix_float* invS);

void ekfNorm(EKF_ctx_t* ctx);
void ekfInitConditions(EKF_ctx_t* ctx, const measures_t* measures);
void qInitEstimate(EKF_ctx_t* ctx, const measures_t* measures);
void PInitEstimate(EKF_ctx_t* ctx);

/**
 * Public Functions Definitions
 */
void ekfInit(EKF_ctx_t* ctx, const measures_t* measures) {
  // Allocate workspace if not already allocated
  if (!ctx->wk) {
    ctx->wk = malloc(sizeof(EKF_work_ctx_t));
  }
  initWorkSpace(ctx);

  ctx->state_current = gsl_vector_float_calloc(STATE_SIZE);
  ctx->state_est = gsl_vector_float_calloc(STATE_SIZE);
  ctx->state_prev = gsl_vector_float_calloc(STATE_SIZE);

  ctx->P_current = gsl_matrix_float_calloc(STATE_SIZE, STATE_SIZE);
  ctx->P_prev = gsl_matrix_float_calloc(STATE_SIZE, STATE_SIZE);
  ctx->P_est = gsl_matrix_float_calloc(STATE_SIZE, STATE_SIZE);

  ctx->acc = gsl_vector_float_calloc(3);
  ctx->velAng = gsl_vector_float_calloc(3);
  ctx->mag = gsl_vector_float_calloc(3);

  ctx->magStdDev = gsl_vector_float_alloc(3);
  gsl_vector_float_set_all(ctx->magStdDev, MAG_STD_DEVIATION);

  ctx->horizonRefG = gsl_vector_float_calloc(3);
  gsl_vector_float_set(ctx->horizonRefG, 0, 0);
  gsl_vector_float_set(ctx->horizonRefG, 1, 0);
  gsl_vector_float_set(ctx->horizonRefG, 2, 1);
  gsl_vector_float_scale(ctx->horizonRefG, GSL_CONST_MKS_GRAV_ACCEL);

  ctx->currentTime = 0;
  ctx->prevTime = 0;

  for (uint8_t i = 0; i < 3; i++) {
    gsl_vector_float_set(ctx->acc, i, measures->acc[i]);
    gsl_vector_float_set(ctx->velAng, i, measures->velAng[i]);
  }

  ekfInitConditions(ctx, measures);
}
void ekfDeinit(EKF_ctx_t* ctx) {
  gsl_quat_float_free(ctx->state_current);
  gsl_quat_float_free(ctx->state_est);
  gsl_quat_float_free(ctx->state_prev);

  gsl_matrix_float_free(ctx->P_current);
  gsl_matrix_float_free(ctx->P_est);
  gsl_matrix_float_free(ctx->P_prev);

  gsl_vector_float_free(ctx->acc);
  gsl_vector_float_free(ctx->velAng);
  gsl_vector_float_free(ctx->mag);

  gsl_vector_float_free(ctx->magStdDev);

  gsl_vector_float_free(ctx->horizonRefG);

  deInitWorkSpace(ctx);
  free(ctx->wk);
}
void ekfStep(EKF_ctx_t* ctx, const measures_t* measures,
             const float currentTime) {
  ekfUpdate(ctx, measures, currentTime);
  ekfEstimate(ctx);
  ekfCorrect(ctx);
  ekfNorm(ctx);
}

/**
 * Private Functions Definitions
 */

void initWorkSpace(EKF_ctx_t* ctx) {
  EKF_work_ctx_t* wk = ctx->wk;

  wk->F = gsl_matrix_float_calloc(STATE_SIZE, STATE_SIZE);
  wk->W = gsl_matrix_float_calloc(STATE_SIZE, CORRECTION_SIZE);
  wk->Q = gsl_matrix_float_calloc(STATE_SIZE, STATE_SIZE);

  wk->H = gsl_matrix_float_calloc(CORRECTION_SIZE, STATE_SIZE);
  wk->K = gsl_matrix_float_calloc(STATE_SIZE, CORRECTION_SIZE);
  wk->R = gsl_matrix_float_calloc(CORRECTION_SIZE, CORRECTION_SIZE);
  wk->S = gsl_matrix_float_calloc(CORRECTION_SIZE, CORRECTION_SIZE);
  wk->invS = gsl_matrix_float_calloc(CORRECTION_SIZE, CORRECTION_SIZE);
  wk->doubleS = gsl_matrix_calloc(CORRECTION_SIZE, CORRECTION_SIZE);
  wk->doubleInvS = gsl_matrix_calloc(CORRECTION_SIZE, CORRECTION_SIZE);

  wk->I4 = gsl_matrix_float_alloc(4, 4);
  gsl_matrix_float_set_identity((gsl_matrix_float*)wk->I4);
  wk->I7 = gsl_matrix_float_alloc(7, 7);
  gsl_matrix_float_set_identity((gsl_matrix_float*)wk->I7);
  wk->I3 = gsl_matrix_float_alloc(3, 3);
  gsl_matrix_float_set_identity((gsl_matrix_float*)wk->I3);

  wk->M1_4_3 = gsl_matrix_float_calloc(4, 3);
  wk->M2_4_4 = gsl_matrix_float_calloc(4, 4);
  wk->M2_7_7 = gsl_matrix_float_calloc(7, 7);
  wk->M1_7_3 = gsl_matrix_float_calloc(7, 3);

  wk->z = gsl_vector_float_calloc(CORRECTION_SIZE);
  wk->h = gsl_vector_float_calloc(CORRECTION_SIZE);

  // NEW: Allocate double-precision buffers for QR inversion (size 6x6)
  wk->inv_tmpMatrix_d = gsl_matrix_alloc(CORRECTION_SIZE, CORRECTION_SIZE);
  wk->inv_tmpTau_d = gsl_vector_alloc(CORRECTION_SIZE);
  wk->inv_tmpB_d = gsl_vector_alloc(CORRECTION_SIZE);
  wk->inv_tmpX_d = gsl_vector_alloc(CORRECTION_SIZE);

  // NEW: Allocate temporary buffers once
  wk->tmpQuat = gsl_quat_float_calloc();
  wk->tmpState = gsl_vector_float_calloc(STATE_SIZE);
  wk->tmp4x4 = gsl_matrix_float_calloc(4, 4);
  wk->tmp7x7 = gsl_matrix_float_calloc(7, 7);
  wk->tmp3x3 = gsl_matrix_float_calloc(3, 3);
  wk->tmpStdDevMat = gsl_matrix_float_calloc(CORRECTION_SIZE, CORRECTION_SIZE);
  wk->tmpBufferMat = gsl_matrix_float_calloc(CORRECTION_SIZE, STATE_SIZE);
  wk->tmp3vec = gsl_vector_float_calloc(3);

  // NEW: Allocate tmpRTransK (6x4) for R*Trans(K)
  wk->tmpRTransK = gsl_matrix_float_calloc(CORRECTION_SIZE, STATE_SIZE);
}

void deInitWorkSpace(EKF_ctx_t* ctx) {
  EKF_work_ctx_t* wk = ctx->wk;

  gsl_matrix_float_free(wk->F);
  gsl_matrix_float_free(wk->W);
  gsl_matrix_float_free(wk->Q);
  gsl_matrix_float_free(wk->K);
  gsl_matrix_float_free(wk->H);
  gsl_matrix_float_free(wk->S);
  gsl_matrix_float_free(wk->invS);
  gsl_matrix_float_free(wk->R);

  gsl_matrix_free(wk->doubleInvS);
  gsl_matrix_free(wk->doubleS);

  gsl_matrix_float_free((gsl_matrix_float*)wk->I7);
  gsl_matrix_float_free((gsl_matrix_float*)wk->I4);
  gsl_matrix_float_free((gsl_matrix_float*)wk->I3);

  gsl_vector_float_free(wk->z);
  gsl_vector_float_free(wk->h);

  gsl_matrix_float_free(wk->M1_4_3);
  gsl_matrix_float_free(wk->M2_4_4);
  gsl_matrix_float_free(wk->M2_7_7);
  gsl_matrix_float_free(wk->M1_7_3);

  // NEW: Free double-precision buffers
  gsl_matrix_free(wk->inv_tmpMatrix_d);
  gsl_vector_free(wk->inv_tmpTau_d);
  gsl_vector_free(wk->inv_tmpB_d);
  gsl_vector_free(wk->inv_tmpX_d);

  // NEW: Free preallocated temporary buffers
  gsl_quat_float_free(wk->tmpQuat);
  gsl_vector_float_free(wk->tmpState);
  gsl_matrix_float_free(wk->tmp4x4);
  gsl_matrix_float_free(wk->tmp7x7);
  gsl_matrix_float_free(wk->tmp3x3);
  gsl_matrix_float_free(wk->tmpStdDevMat);
  gsl_matrix_float_free(wk->tmpBufferMat);
  gsl_vector_float_free(wk->tmp3vec);

  // NEW: Free tmpRTransK
  gsl_matrix_float_free(wk->tmpRTransK);
}

void ekfUpdate(EKF_ctx_t* ctx, const measures_t* measures,
               const float currentTime) {
  ctx->prevTime = ctx->currentTime;
  ctx->currentTime = currentTime;

  gsl_vector_float_memcpy(ctx->state_prev, ctx->state_current);
  gsl_matrix_float_memcpy(ctx->P_prev, ctx->P_current);

  for (uint8_t i = 0; i < 3; i++) {
    gsl_vector_float_set(ctx->acc, i, measures->acc[i]);
    gsl_vector_float_set(ctx->velAng, i, measures->velAng[i]);
    gsl_vector_float_set(ctx->mag, i, measures->mag[i]);
  }

  // Normalize acc
  // float accNorm = 0;
  // gsl_blas_sdot(ctx->acc, ctx->acc, &accNorm);
  // accNorm = sqrtf(accNorm);
  // gsl_vector_float_scale(ctx->acc, ACC_SCALE / accNorm);

  // Normalize mag
  float magNorm = 0;
  gsl_blas_sdot(ctx->mag, ctx->mag, &magNorm);
  magNorm = sqrtf(magNorm);
  gsl_vector_float_scale(ctx->mag, MAG_SCALE / magNorm);
}

void ekfEstimate(EKF_ctx_t* ctx) {
  qEst(ctx);
  PEst(ctx);
}

// Modified qEst: use tmpQuat and tmp4x4 instead of per-call allocations.
void qEst(EKF_ctx_t* ctx) {
  EKF_work_ctx_t* wk = ctx->wk;
  float deltaT = ctx->currentTime - ctx->prevTime;
  qEstPrimitive(ctx->velAng, deltaT, ctx->state_prev, ctx->state_est,
                wk->tmpQuat, wk->tmp4x4, wk->M2_7_7);
  gsl_vector_float_view quat_est = gsl_vector_float_subvector(ctx->state_est,0,4);
  gsl_quat_float_normalize(&quat_est.vector);
}

void PEst(EKF_ctx_t* ctx) {
  EKF_work_ctx_t* wk = ctx->wk;
  getF(ctx);
  getQ(ctx);
  // Here W is not used in the primitive implementation.
  PEstPrimitive(ctx->P_prev, wk->F, wk->Q, ctx->P_est, wk->tmp7x7);
  
}

// Modify getF to use tmpQuat (from workspace) instead of allocating a new one.
void getF(EKF_ctx_t* ctx) {
  EKF_work_ctx_t* wk = ctx->wk;
  // Use workspace tmpQuat instead of a new allocation.
  getFPrimitive(ctx->velAng, ctx->currentTime - ctx->prevTime, ctx->state_prev,
                wk->tmpQuat, ctx->wk->F, wk->I7);
}

void getW(EKF_ctx_t* ctx) {
  getWPrimitive(ctx->state_prev, ctx->currentTime - ctx->prevTime, ctx->wk->W);
}

// Similarly update getQ, get_h, and getH to use wk->tmpStdDevMat,
// wk->tmpBufferMat, wk->tmp3x3, and wk->tmp3vec. For example, in getQ:
void getQ(EKF_ctx_t* ctx) {
  getW(ctx);

  gsl_blas_sgemm(CblasNoTrans, CblasTrans, GYRO_STD_DEVIATION, ctx->wk->W,
                 ctx->wk->W, 0, ctx->wk->Q);
  gsl_matrix_float_set(ctx->wk->Q, 4, 4, GYRO_BIAS_NOISE);
  gsl_matrix_float_set(ctx->wk->Q, 5, 5, GYRO_BIAS_NOISE);
  gsl_matrix_float_set(ctx->wk->Q, 6, 6, GYRO_BIAS_NOISE);
}

void ekfCorrect(EKF_ctx_t* ctx) {
  EKF_work_ctx_t* wk = ctx->wk;
  gsl_vector_float* z = wk->z;
  gsl_vector_float* h = wk->h;
  for (uint8_t i = 0; i < 3; i++) {
    gsl_vector_float_set(z, i, gsl_vector_float_get(ctx->acc, i));
  }
  get_h(ctx);
  gsl_vector_float_sub(z, h);

  getK(ctx);
  // Replace the in‐line update by a call to the primitive.
  qCorrectPrimitive(ctx->state_est, wk->K, z, ctx->state_current, wk->tmpState);
  PCorrect(ctx);

  if (ctx->currentTime - ctx->lastMagCorrectionTime_s >
      MAG_CORRECTION_PERIOD_S) {
    ctx->lastMagCorrectionTime_s = ctx->currentTime;
    correctMag(ctx->P_current, ctx->state_current, ctx->mag, ctx->magStdDev);
  }
}

void PCorrect(EKF_ctx_t* ctx) {
  EKF_work_ctx_t* wk = ctx->wk;
  PCorrectPrimitive(ctx->P_est, wk->K, wk->H, wk->R, ctx->P_current, wk->tmp7x7,
                    wk->tmpRTransK, wk->I7);
}

void get_h(EKF_ctx_t* ctx) {
  EKF_work_ctx_t* wk = ctx->wk;
  get_hPrimitive(ctx->state_est, ctx->horizonRefG, wk->h);
}

int getH(EKF_ctx_t* ctx) {
  EKF_work_ctx_t* wk = ctx->wk;
  gsl_matrix_float_set_zero(wk->H);
  // Use two columns from tmpBufferMat as temporary vectors for pV1 and pV2.
  gsl_vector_float_view tmp_v1 = gsl_matrix_float_column(wk->tmpBufferMat, 0);
  gsl_vector_float_view tmp_v2 = gsl_matrix_float_column(wk->tmpBufferMat, 1);
  gsl_vector_float_view quat_est =
      gsl_vector_float_subvector(ctx->state_est, 0, 4);
  // Delegate to the primitive using workspace buffers:
  getHPrimitive((gsl_quat_float*)&quat_est.vector, ctx->horizonRefG, wk->H,
                wk->tmp3vec, wk->tmp3x3, &tmp_v1.vector, &tmp_v2.vector);
  return 0;
}

int getR(EKF_ctx_t* ctx) {
  getRPrimitive(ctx->wk->R);
  return 0;
}

// S = H * P_est * trans(H) + R
int getS(EKF_ctx_t* ctx) {
  getR(ctx);
  getSPrimitive(ctx->wk->H, ctx->P_est, ctx->wk->R, ctx->wk->S,
                ctx->wk->M1_7_3);
  return 0;
}

// K = P_est * trans(H) * inv(S)
int getK(EKF_ctx_t* ctx) {
  getH(ctx);
  getS(ctx);
  invertMatrixFloat(ctx, ctx->wk->S, ctx->wk->invS);
  getKPrimitive(ctx->P_est, ctx->wk->H, ctx->wk->invS, ctx->wk->K,
                ctx->wk->M1_7_3);
  return 0;
}

void invertMatrixFloatPrimitive(const gsl_matrix* S, gsl_matrix* invS,
                                gsl_vector* tau, gsl_vector* b, gsl_vector* x) {
  int n = S->size1;

  // Create a temporary copy of S
  gsl_matrix* temp = gsl_matrix_alloc(S->size1, S->size2);
  gsl_matrix_memcpy(temp, S);

  // Perform QR decomposition on the temporary copy
  gsl_linalg_QR_decomp(temp, tau);

  // Solve for each column of the identity matrix
  for (int i = 0; i < n; i++) {
    gsl_vector_set_zero(b);
    gsl_vector_set(b, i, 1.0);
    gsl_linalg_QR_solve(temp, tau, b, x);

    // Copy the solution into the inverse matrix
    for (int j = 0; j < n; j++) {
      gsl_matrix_set(invS, j, i, gsl_vector_get(x, j));
    }
  }

  // Free the temporary matrix
  gsl_matrix_free(temp);
}

void invertMatrixFloat(EKF_ctx_t* ctx, const gsl_matrix_float* S,
                       gsl_matrix_float* invS) {
  int n = S->size1;
  EKF_work_ctx_t* wk = ctx->wk;

  // Convert float matrix S into double-precision matrix using preallocated
  // buffer
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      gsl_matrix_set(wk->inv_tmpMatrix_d, i, j, gsl_matrix_float_get(S, i, j));
    }
  }

  // Call the primitive function to compute the inverse
  invertMatrixFloatPrimitive(wk->inv_tmpMatrix_d, wk->doubleInvS,
                             wk->inv_tmpTau_d, wk->inv_tmpB_d, wk->inv_tmpX_d);

  // Copy the double solution into the float output matrix
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      gsl_matrix_float_set(invS, j, i,
                           (float)gsl_matrix_get(wk->doubleInvS, j, i));
    }
  }
}

void ekfNorm(EKF_ctx_t* ctx) {
  gsl_vector_float_view q_current =
      gsl_vector_float_subvector(ctx->state_current, 0, 4);
  gsl_quat_float_normalize(&q_current.vector);
}

void ekfInitConditions(EKF_ctx_t* ctx, const measures_t* measures) {
  qInitEstimate(ctx, measures);
  PInitEstimate(ctx);
}

void qInitEstimate(EKF_ctx_t* ctx, const measures_t* measures) {
  gsl_quat_float_set(ctx->state_current, 0, 1);
  gsl_quat_float_set(ctx->state_current, 1, 0);
  gsl_quat_float_set(ctx->state_current, 2, 0);
  gsl_quat_float_set(ctx->state_current, 3, 0);
}

void PInitEstimate(EKF_ctx_t* ctx) {
  gsl_matrix_float_set_identity(ctx->P_current);
}

// ...existing code...

void qEstPrimitive(const gsl_vector_float* velAng, float deltaT,
                   const gsl_quat_float* statePrev, gsl_quat_float* stateEst,
                   gsl_quat_float* tmpQuat, gsl_matrix_float* qVelAngMat,
                   gsl_matrix_float* state1Mat) {
  // Removed dynamic allocation; use provided tmpQuat, qVelAngMat, and q1Mat
  // instead.
  gsl_quat_float_set(tmpQuat, 0, 0);
  gsl_quat_float_set(
      tmpQuat, 1,
      gsl_vector_float_get(velAng, 0) + gsl_vector_float_get(statePrev, 4));

  gsl_quat_float_set(
      tmpQuat, 2,
      gsl_vector_float_get(velAng, 1) + gsl_vector_float_get(statePrev, 5));
  gsl_quat_float_set(
      tmpQuat, 3,
      gsl_vector_float_get(velAng, 2) + gsl_vector_float_get(statePrev, 6));

  gsl_quat_float_toMatrix(tmpQuat, qVelAngMat);
  gsl_matrix_float_set_identity(state1Mat);
  gsl_matrix_float_scale(qVelAngMat, deltaT / 2.f);

  gsl_matrix_float_view q1Mat =
      gsl_matrix_float_submatrix(state1Mat, 0, 0, 4, 4);

  gsl_matrix_float_add(&q1Mat.matrix, qVelAngMat);
  gsl_blas_sgemv(CblasNoTrans, 1.F, state1Mat, statePrev, 0, stateEst);
}

void PEstPrimitive(const gsl_matrix_float* PPrev, const gsl_matrix_float* F,
                   const gsl_matrix_float* Q, gsl_matrix_float* PEst,
                   gsl_matrix_float* tmp7x7) {
  // Use provided tmp4x4 instead of allocating a new one.
  gsl_blas_sgemm(CblasNoTrans, CblasTrans, 1, PPrev, F, 0, tmp7x7);
  gsl_blas_sgemm(CblasNoTrans, CblasNoTrans, 1, F, tmp7x7, 0, PEst);
  gsl_matrix_float_scale(PEst, P_ESTIMATE_SCALE);
  gsl_matrix_float_add(PEst, Q);
}

void PCorrectPrimitive(const gsl_matrix_float* P, const gsl_matrix_float* K,
                       const gsl_matrix_float* H, const gsl_matrix_float* R,
                       gsl_matrix_float* PCorrect, gsl_matrix_float* tmp4x4,
                       gsl_matrix_float* tmp6x4, gsl_matrix_float* I4) {
  // Set I4 to identity (caller may preallocate and reuse I4)
  gsl_matrix_float_set_identity(I4);
  // Compute I - K*H into tmp4x4
  gsl_blas_sgemm(CblasNoTrans, CblasNoTrans, -1, K, H, 0, tmp4x4);
  gsl_matrix_float_add(tmp4x4, I4);
  gsl_matrix_float* tmpPCorrect = gsl_matrix_float_calloc(P->size1, P->size1);

  if (P_CORRECT_METHOD == 0) {
    gsl_blas_sgemm(CblasNoTrans, CblasNoTrans, 1, tmp4x4, P, 0, PCorrect);

  } else if (P_CORRECT_METHOD == 1) {
    gsl_blas_sgemm(CblasNoTrans, CblasTrans, 1, P, tmp4x4, 0, tmpPCorrect);
    gsl_blas_sgemm(CblasNoTrans, CblasNoTrans, 1, tmp4x4, tmpPCorrect, 0,
                   PCorrect);

    gsl_blas_sgemm(CblasNoTrans, CblasTrans, 1, R, K, 0, tmp6x4);
    gsl_blas_sgemm(CblasNoTrans, CblasNoTrans, 1, K, tmp6x4, 0, tmp4x4);
    gsl_matrix_float_add(PCorrect, tmp4x4);
  }
  gsl_matrix_float_free(tmpPCorrect);
}

void qCorrectPrimitive(const gsl_quat_float* q, const gsl_matrix_float* K,
                       const gsl_vector_float* h, gsl_quat_float* qCorrect,
                       gsl_vector_float* tmp) {
  // Use the provided tmp instead of dynamic allocation.
  gsl_blas_sgemv(CblasNoTrans, 1, K, h, 0, tmp);  // tmp = K*h
  gsl_vector_float_memcpy(qCorrect, q);
  gsl_vector_float_add(qCorrect, tmp);  // qCorrect = q + K*h
}

void getKPrimitive(const gsl_matrix_float* P, const gsl_matrix_float* H,
                   const gsl_matrix_float* invS, gsl_matrix_float* K,
                   gsl_matrix_float* tmp4x6) {
  // Use provided tmp4x6 instead of allocating it.
  gsl_blas_sgemm(CblasTrans, CblasNoTrans, 1, H, invS, 0,
                 tmp4x6);  // tmp4x6 = trans(H)*invS
  gsl_blas_sgemm(CblasNoTrans, CblasNoTrans, 1, P, tmp4x6, 0,
                 K);  // K = P * tmp4x6
}

void getHPrimitive(const gsl_vector_float* state_est, const gsl_vector_float* acc,
                   gsl_matrix_float* H, gsl_vector_float* pQv,
                   gsl_matrix_float* pM2, gsl_vector_float* pV1,
                   gsl_vector_float* pV2) {
  gsl_matrix_float_set_zero(H);

  gsl_vector_float_view q_view = gsl_vector_float_subvector(state_est,0,4);
  gsl_vector_float *q = &q_view.vector;

  gsl_quat_float_get_imaginary(q, pQv);

  // Upper part (acc related)
  gsl_vector_float_view u_g = gsl_matrix_float_subcolumn(H, 0, 0, 3);
  cross_product(acc, pQv, &u_g.vector);

  gsl_matrix_float_view upperRightH = gsl_matrix_float_submatrix(H, 0, 1, 3, 3);
  gsl_matrix_float* pM1 = &upperRightH.matrix;

  gsl_vector_float_memcpy(pV2, acc);
  gsl_vector_float_scale(pV2, gsl_quat_float_get(q, 0));  // scale by qw
  gsl_vector_float_memcpy(pV1, &u_g.vector);
  gsl_vector_float_add(pV1, pV2);  // ug + qw*acc

  skewSymFromVector(pV1, pM1);  // [ug + qw*acc]x

  gsl_matrix_float_set_zero(pM2);

  gsl_blas_sger(1, acc, pQv, pM2);  // pM2 = acc ⊗ pQv
  gsl_matrix_float_scale(pM2, -1.f);
  gsl_matrix_float_add(pM1, pM2);  // [ug + qw*acc]x - acc⊗pQv

  gsl_matrix_float_set_identity(pM2);
  float acc_dot = 0;
  gsl_blas_sdot(pQv, acc, &acc_dot);
  gsl_matrix_float_scale(pM2, acc_dot);
  gsl_matrix_float_add(pM1, pM2);  // finalize upper part

  gsl_matrix_float_scale(H, 2.f);
}

void getRPrimitive(gsl_matrix_float* R) {
  gsl_matrix_float_set_zero(R);
  for (int i = 0; i < 3; i++) {
    gsl_matrix_float_set(R, i, i, ACC_STD_DEVIATION);
  }
}

void getSPrimitive(const gsl_matrix_float* H, const gsl_matrix_float* P,
                   const gsl_matrix_float* R, gsl_matrix_float* S,
                   gsl_matrix_float* tmp4x6) {
  gsl_matrix_float_set_zero(S);
  // Use provided tmp4x6 instead of allocating it.
  gsl_blas_sgemm(CblasNoTrans, CblasTrans, 1, P, H, 0,
                 tmp4x6);  // tmp4x6 = P*trans(H)
  gsl_blas_sgemm(CblasNoTrans, CblasNoTrans, 1, H, tmp4x6, 0,
                 S);           // S = H*tmp4x6
  gsl_matrix_float_add(S, R);  // S = H*P*trans(H) + R
}

void getFPrimitive(const gsl_vector_float* velAng, float deltaT,
                   gsl_vector_float* prevState, gsl_quat_float* tmpQuat,
                   gsl_matrix_float* F, const gsl_matrix_float* I7) {
  gsl_matrix_float_set_zero(F);
  // Use workspace tmpQuat instead of a new allocation.
  gsl_vector_float_view estQuat = gsl_vector_float_subvector(prevState, 0, 4);
  gsl_vector_float_view estBias = gsl_vector_float_subvector(prevState, 4, 3);
  gsl_matrix_float_view quatsF = gsl_matrix_float_submatrix(F, 0, 0, 4, 4);

  gsl_vector_float* velAngNoBias = gsl_vector_float_alloc(3);
  gsl_vector_float_memcpy(velAngNoBias, velAng);
  gsl_vector_float_add(velAngNoBias, &estBias.vector);

  gsl_quat_float_fromVector(velAngNoBias, tmpQuat);
  gsl_quat_float_toMatrix(tmpQuat, &quatsF.matrix);

  float q1, q2, q3, q4;
  q1 = gsl_quat_float_get(&estQuat.vector, 0);
  q2 = gsl_quat_float_get(&estQuat.vector, 1);
  q3 = gsl_quat_float_get(&estQuat.vector, 2);
  q4 = gsl_quat_float_get(&estQuat.vector, 3);

  gsl_matrix_float_set(F, 0, 4, -q2);
  gsl_matrix_float_set(F, 0, 5, -q3);
  gsl_matrix_float_set(F, 0, 6, -q4);

  gsl_matrix_float_set(F, 1, 4, q1);
  gsl_matrix_float_set(F, 1, 5, -q4);
  gsl_matrix_float_set(F, 1, 6, -q3);

  gsl_matrix_float_set(F, 2, 4, q4);
  gsl_matrix_float_set(F, 2, 5, q1);
  gsl_matrix_float_set(F, 2, 6, -q2);

  gsl_matrix_float_set(F, 3, 4, -q3);
  gsl_matrix_float_set(F, 3, 5, q2);
  gsl_matrix_float_set(F, 3, 6, q1);

  gsl_matrix_float_scale(F, (deltaT) * 0.5f);
  gsl_matrix_float_add(F, I7);

  gsl_vector_float_free(velAngNoBias);
}

void getWPrimitive(const gsl_quat_float* qPrev, float deltaT,
                   gsl_matrix_float* W) {
  float q1, q2, q3, q4;
  q1 = gsl_quat_float_get(qPrev, 0);
  q2 = gsl_quat_float_get(qPrev, 1);
  q3 = gsl_quat_float_get(qPrev, 2);
  q4 = gsl_quat_float_get(qPrev, 3);
  gsl_matrix_float_set_zero(W);

  gsl_matrix_float_set(W, 0, 0, -q2);
  gsl_matrix_float_set(W, 0, 1, -q3);
  gsl_matrix_float_set(W, 0, 2, -q4);

  gsl_matrix_float_set(W, 1, 0, q1);
  gsl_matrix_float_set(W, 1, 1, -q4);
  gsl_matrix_float_set(W, 1, 2, q3);

  gsl_matrix_float_set(W, 2, 0, q4);
  gsl_matrix_float_set(W, 2, 1, q1);
  gsl_matrix_float_set(W, 2, 2, -q2);

  gsl_matrix_float_set(W, 3, 0, -q3);
  gsl_matrix_float_set(W, 3, 1, q2);
  gsl_matrix_float_set(W, 3, 2, q1);

  gsl_matrix_float_scale(W, (deltaT) * 0.5f);
}

void get_hPrimitive(const gsl_vector_float* state_est,
                    const gsl_vector_float* horizonRefG, gsl_vector_float* h) {
  gsl_matrix_float* pRotMat = gsl_matrix_float_calloc(3, 3);
  gsl_vector_float_view ExpectedAcc = gsl_vector_float_subvector(h, 0, 3);
  gsl_vector_float_view q_est = gsl_vector_float_subvector(state_est, 0, 4);

  gsl_quat_float_toRotMatrix(&q_est.vector, pRotMat);

  gsl_blas_sgemv(CblasTrans, 1, pRotMat, horizonRefG, 0, &ExpectedAcc.vector);

  gsl_matrix_float_free(pRotMat);
}