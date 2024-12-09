#include "EKF.h"

#include "EKF_const.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_math.h"
#include "gsl/gsl_matrix_float.h"
#include "gsl/gsl_vector_float.h"
#include "gsl_quaternion_float.h"
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

void initWorkSpace(EKF_ctx_t *ctx);
void deInitWorkSpace(EKF_ctx_t *ctx);
void ekfUpdate(EKF_ctx_t *ctx, const measures_t *measures,
               const float currentTime);

void ekfEstimate(EKF_ctx_t *ctx);
void qEst(EKF_ctx_t *ctx);
void PEst(EKF_ctx_t *ctx);
void getF(EKF_ctx_t *ctx);
void getQ(EKF_ctx_t *ctx);
void getW(EKF_ctx_t *ctx);

void ekfCorrect(EKF_ctx_t *ctx);
void get_h(EKF_ctx_t *ctx);
int getH(EKF_ctx_t *ctx);
int getR(EKF_ctx_t *ctx);
int getS(EKF_ctx_t *ctx);
int getK(EKF_ctx_t *ctx);

void ekfNorm(EKF_ctx_t *ctx);
void ekfInitConditions(EKF_ctx_t *ctx, const measures_t *measures);
void qInitEstimate(EKF_ctx_t *ctx, const measures_t *measures);
void PInitEstimate(EKF_ctx_t *ctx);

/**
 * Public Functions Definitions
 */
void ekfInit(EKF_ctx_t *ctx, EKF_work_ctx_t *wk_ctx,
             const measures_t *measures) {
  ctx->wk = wk_ctx;
  ctx->q_current = gsl_quat_float_calloc();
  ctx->q_est = gsl_quat_float_calloc();
  ctx->q_prev = gsl_quat_float_calloc();

  ctx->P_current = gsl_matrix_float_calloc(P_SIZE, P_SIZE);
  ctx->P_prev = gsl_matrix_float_calloc(P_SIZE, P_SIZE);
  ctx->P_est = gsl_matrix_float_calloc(P_SIZE, P_SIZE);

  ctx->acc = gsl_vector_float_calloc(3);
  ctx->mag = gsl_vector_float_calloc(3);
  ctx->velAng = gsl_vector_float_calloc(3);

  ctx->latitude = LATITUDE_DEG * M_PI / 180;

  ctx->horizonRefG = gsl_vector_float_calloc(3);
  gsl_vector_float_set(ctx->horizonRefG, 0, 0);
  gsl_vector_float_set(ctx->horizonRefG, 1, 0);
  gsl_vector_float_set(ctx->horizonRefG, 2, 9.8);

  ctx->horizonRefMag = gsl_vector_float_calloc(3);
  gsl_vector_float_set(ctx->horizonRefMag, 0, cosf(ctx->latitude));
  gsl_vector_float_set(ctx->horizonRefMag, 1, 0);
  gsl_vector_float_set(ctx->horizonRefMag, 2, -sinf(ctx->latitude));

  for (uint8_t i = 0; i < 3; i++) {
    gsl_vector_float_set(ctx->acc, i, measures->acc[i]);
    gsl_vector_float_set(ctx->mag, i, measures->mag[i]);
    gsl_vector_float_set(ctx->velAng, i, measures->velAng[i]);
  }

  initWorkSpace(ctx);
  ekfInitConditions(ctx, measures);
}
void ekfDeinit(EKF_ctx_t *ctx) {
  gsl_quat_float_free(ctx->q_current);
  gsl_quat_float_free(ctx->q_est);
  gsl_quat_float_free(ctx->q_prev);

  gsl_matrix_float_free(ctx->P_current);
  gsl_matrix_float_free(ctx->P_est);
  gsl_matrix_float_free(ctx->P_prev);

  gsl_quat_float_free(ctx->acc);
  gsl_quat_float_free(ctx->mag);
  gsl_quat_float_free(ctx->velAng);

  gsl_quat_float_free(ctx->horizonRefG);
  gsl_quat_float_free(ctx->horizonRefMag);

  deInitWorkSpace(ctx);
}
void ekfStep(EKF_ctx_t *ctx, const measures_t *measures,
             const float currentTime) {
  ekfUpdate(ctx, measures, currentTime);
  ekfEstimate(ctx);
  ekfCorrect(ctx);
  ekfNorm(ctx);
}

/**
 * Private Functions Definitions
 */

void initWorkSpace(EKF_ctx_t *ctx) {
  EKF_work_ctx_t *wk = ctx->wk;

  wk->F = gsl_matrix_float_calloc(4, 4);
  wk->W = gsl_matrix_float_calloc(4, 3);
  wk->Q = gsl_matrix_float_calloc(4, 4);

  wk->H = gsl_matrix_float_calloc(6, 4);
  wk->K = gsl_matrix_float_calloc(4, 6);
  wk->R = gsl_matrix_float_calloc(6, 6);
  wk->S = gsl_matrix_float_calloc(6, 6);
  wk->invS = gsl_matrix_float_calloc(6, 6);
  wk->doubleS = gsl_matrix_calloc(6, 6);
  wk->doubleInvS = gsl_matrix_calloc(6, 6);

  wk->I4 = gsl_matrix_float_alloc(4, 4);
  gsl_matrix_float_set_identity((gsl_matrix_float *)wk->I4);
  wk->I3 = gsl_matrix_float_alloc(3, 3);
  gsl_matrix_float_set_identity((gsl_matrix_float *)wk->I3);

  wk->M1_4_6 = gsl_matrix_float_calloc(4, 6);
  wk->M2_4_4 = gsl_matrix_float_calloc(4, 4);

  wk->z = gsl_vector_float_calloc(6);
  wk->h = gsl_vector_float_calloc(6);
}

void deInitWorkSpace(EKF_ctx_t *ctx) {
  EKF_work_ctx_t *wk = ctx->wk;

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

  gsl_matrix_float_free((gsl_matrix_float *)wk->I4);
  gsl_matrix_float_free((gsl_matrix_float *)wk->I3);

  gsl_vector_float_free(wk->z);
  gsl_vector_float_free(wk->h);

  // gsl_matrix_float_free(wk->M1_4_6);
  // gsl_matrix_float_free(wk->M2_4_4);
  // gsl_vector_float_free(wk->v1);
  // gsl_vector_float_free(wk->v2);
}

void ekfUpdate(EKF_ctx_t *ctx, const measures_t *measures,
               const float currentTime) {
  ctx->prevTime = ctx->currentTime;
  ctx->currentTime = currentTime;

  gsl_vector_float_memcpy(ctx->q_prev, ctx->q_current);
  gsl_vector_float_set_zero(ctx->q_current);
  gsl_vector_float_set_zero(ctx->q_est);

  gsl_matrix_float_memcpy(ctx->P_prev, ctx->P_current);
  gsl_matrix_float_set_zero(ctx->P_current);
  gsl_matrix_float_set_zero(ctx->P_est);

  for (uint8_t i = 0; i < 3; i++) {
    gsl_vector_float_set(ctx->acc, i, measures->acc[i]);
    gsl_vector_float_set(ctx->mag, i, measures->mag[i]);
    gsl_vector_float_set(ctx->velAng, i, measures->velAng[i]);
  }
}

void ekfEstimate(EKF_ctx_t *ctx) {
  qEst(ctx);
  PEst(ctx);
}
void qEst(EKF_ctx_t *ctx) {
  gsl_quat_float *qVelAng = gsl_quat_float_alloc();
  gsl_matrix_float *q1Mat = gsl_matrix_float_alloc(4, 4);

  gsl_quat_float_set(qVelAng, 0, 0);
  gsl_quat_float_set(qVelAng, 1, gsl_vector_float_get(ctx->velAng, 0));
  gsl_quat_float_set(qVelAng, 2, gsl_vector_float_get(ctx->velAng, 1));
  gsl_quat_float_set(qVelAng, 3, gsl_vector_float_get(ctx->velAng, 2));

  gsl_matrix_float *qVelAngMat = gsl_matrix_float_alloc(4, 4);
  gsl_quat_float_toMatrix(qVelAng, qVelAngMat);

  gsl_matrix_float_set_identity(q1Mat);

  gsl_matrix_float_scale(qVelAngMat, (ctx->currentTime - ctx->prevTime) / 2.f);
  gsl_matrix_float_add(q1Mat, qVelAngMat);

  gsl_blas_sgemv(CblasNoTrans, 1.F, q1Mat, ctx->q_prev, 0, ctx->q_est);

  gsl_matrix_float_free(qVelAngMat);
  gsl_matrix_float_free(q1Mat);
  gsl_quat_float_free(qVelAng);
}

void PEst(EKF_ctx_t *ctx) {
  gsl_matrix_float *F = ctx->wk->F;
  gsl_matrix_float *Q = ctx->wk->Q;
  gsl_matrix_float *M2_4_4 = ctx->wk->M2_4_4;

  getF(ctx);
  getQ(ctx);

  gsl_blas_sgemm(CblasNoTrans, CblasTrans, 1, ctx->P_prev, F, 0, M2_4_4);
  gsl_blas_sgemm(CblasNoTrans, CblasNoTrans, 1, F, M2_4_4, 0, ctx->P_est);
  
  gsl_matrix_float_add(ctx->P_est, Q);
}

void getF(EKF_ctx_t *ctx) {
  gsl_matrix_float *F = ctx->wk->F;
  gsl_quat_float *qVelAng = gsl_quat_float_alloc();
  gsl_quat_float_fromVector(ctx->velAng, qVelAng);
  gsl_quat_float_toMatrix(qVelAng, F);

  const gsl_matrix_float *pI = ctx->wk->I4;
  gsl_matrix_float_scale(F, (ctx->currentTime - ctx->prevTime) / 2);
  gsl_matrix_float_add(F, pI);

  gsl_quat_float_free(qVelAng);
}
void getW(EKF_ctx_t *ctx) {
  gsl_matrix_float *W = ctx->wk->W;
  gsl_matrix_float_set_zero(W);

  gsl_matrix_float_set(W, 0, 0, -gsl_quat_float_get(ctx->q_prev, 1));
  gsl_matrix_float_set(W, 0, 1, -gsl_quat_float_get(ctx->q_prev, 2));
  gsl_matrix_float_set(W, 0, 2, -gsl_quat_float_get(ctx->q_prev, 3));

  gsl_matrix_float_set(W, 1, 0, gsl_quat_float_get(ctx->q_prev, 0));
  gsl_matrix_float_set(W, 1, 1, -gsl_quat_float_get(ctx->q_prev, 3));
  gsl_matrix_float_set(W, 1, 2, gsl_quat_float_get(ctx->q_prev, 2));

  gsl_matrix_float_set(W, 2, 0, gsl_quat_float_get(ctx->q_prev, 3));
  gsl_matrix_float_set(W, 2, 1, gsl_quat_float_get(ctx->q_prev, 0));
  gsl_matrix_float_set(W, 2, 2, -gsl_quat_float_get(ctx->q_prev, 1));

  gsl_matrix_float_set(W, 3, 0, -gsl_quat_float_get(ctx->q_prev, 2));
  gsl_matrix_float_set(W, 3, 1, gsl_quat_float_get(ctx->q_prev, 1));
  gsl_matrix_float_set(W, 3, 2, gsl_quat_float_get(ctx->q_prev, 0));

  gsl_matrix_float_scale(W, (ctx->currentTime - ctx->prevTime) / 2.f);
}
void getQ(EKF_ctx_t *ctx) {
  gsl_matrix_float *pStdDevMat = gsl_matrix_float_alloc(3, 3);
  gsl_matrix_float *pBufferMat = gsl_matrix_float_calloc(3, 4);
  getW(ctx);
  gsl_matrix_float_set_identity(pStdDevMat);
  gsl_matrix_float_scale(pStdDevMat, GYRO_STD_DEVIATION);
  gsl_blas_sgemm(CblasNoTrans, CblasTrans, 1, pStdDevMat, ctx->wk->W, 0,
                 pBufferMat);
  gsl_blas_sgemm(CblasNoTrans, CblasNoTrans, 1, ctx->wk->W, pBufferMat, 0,
                 ctx->wk->Q);

  gsl_matrix_float_free(pStdDevMat);
  gsl_matrix_float_free(pBufferMat);
}

void ekfCorrect(EKF_ctx_t *ctx) {
  gsl_vector_float *z = ctx->wk->z;
  gsl_vector_float *v;
  for (uint8_t i = 0; i < 3; i++) {
    gsl_vector_float_set(z, i, gsl_vector_float_get(ctx->acc, i));
    gsl_vector_float_set(z, i + 3, gsl_vector_float_get(ctx->mag, i));
  }

  gsl_vector_float *h = ctx->wk->h;
  get_h(ctx);
  gsl_vector_float_sub(z, h);
  v = z;
  gsl_matrix_float *K = ctx->wk->K;
  getK(ctx);

  gsl_blas_sgemv(CblasNoTrans, 1, K, v, 0, ctx->q_current);
  gsl_vector_float_add(ctx->q_current, ctx->q_est);

  gsl_blas_sgemm(CblasNoTrans, CblasNoTrans, -1, K, ctx->wk->H, 0, ctx->P_prev);

  gsl_matrix_float_add(ctx->P_prev, ctx->wk->I4);

  gsl_blas_sgemm(CblasNoTrans, CblasNoTrans, 1, ctx->P_prev, ctx->P_est, 0,
                 ctx->P_current);
}

void get_h(EKF_ctx_t *ctx) {
  gsl_vector_float *h = ctx->wk->h;
  gsl_matrix_float *pRotMat = gsl_matrix_float_calloc(3, 3);
  gsl_vector_float_view ExpectedAcc = gsl_vector_float_subvector(h, 0, 3);
  gsl_vector_float_view ExpectedMag = gsl_vector_float_subvector(h, 3, 3);

  gsl_quat_float_toRotMatrix(ctx->q_est, pRotMat);

  gsl_blas_sgemv(CblasTrans, 1, pRotMat, ctx->horizonRefG, 0,
                 &ExpectedAcc.vector);
  gsl_blas_sgemv(CblasTrans, 1, pRotMat, ctx->horizonRefMag, 0,
                 &ExpectedMag.vector);

  gsl_matrix_float_free(pRotMat);
}

/*
//[ug+qw*r, skewSym(ur+qw*r)+dot(qv,r).*eye(3)]
int getH(EKF_ctx_t *ctx) {
  gsl_matrix_float *H = ctx->wk->H;
  gsl_matrix_float_set_zero(H);
  
  gsl_vector_float_view u_g = gsl_matrix_float_subcolumn(H, 0, 0, 3);
  gsl_vector_float_view u_r = gsl_matrix_float_subcolumn(H, 0, 3, 3);


  gsl_matrix_float_view upperRightH = gsl_matrix_float_submatrix(H, 0, 1, 3, 3);
  gsl_matrix_float_view lowerRightH = gsl_matrix_float_submatrix(H, 3, 1, 3, 3);

  gsl_vector_float *pQv = gsl_vector_float_alloc(3);
  gsl_quat_float_get_imaginary(ctx->q_est, pQv);

  cross_product(ctx->horizonRefG, pQv, &u_g.vector);
  cross_product(ctx->horizonRefMag, pQv, &u_r.vector);

  gsl_vector_float *pV1 = gsl_vector_float_alloc(3);

  gsl_vector_float_memcpy(pV1, ctx->horizonRefG);
  gsl_vector_float_scale(pV1, gsl_quat_float_get(ctx->q_est, 0));
  gsl_vector_float_add(&u_g.vector, pV1);

  gsl_vector_float_memcpy(pV1, ctx->horizonRefMag);
  gsl_vector_float_scale(pV1, gsl_quat_float_get(ctx->q_est, 0));
  gsl_vector_float_add(&u_r.vector, pV1);

  // upperRightH

  gsl_matrix_float *pM1 = &upperRightH.matrix;
  gsl_matrix_float *pM2 = gsl_matrix_float_calloc(3, 3);
  gsl_vector_float *pV2 = gsl_vector_float_alloc(3);
  


  skewSymFromVector(&u_g.vector, pM1);

  gsl_matrix_float_set_identity(pM2);

  float G_Qv_dot = 0;
  gsl_blas_sdot(pQv, ctx->horizonRefG, &G_Qv_dot);
  gsl_matrix_float_scale(pM2, G_Qv_dot);

  gsl_matrix_float_add(pM1, pM2);

  // lowerRightH

  pM1 = &lowerRightH.matrix;
  gsl_matrix_float_set_zero(pM2);

  skewSymFromVector(&u_r.vector, pM1);


  gsl_matrix_float_set_identity(pM2);

  float Mag_Qv_dot = 0;
  gsl_blas_sdot(pQv, ctx->horizonRefMag, &Mag_Qv_dot);
  gsl_matrix_float_scale(pM2, Mag_Qv_dot);

  gsl_matrix_float_add(pM1, pM2);
  gsl_matrix_float_scale(H, 2.f);

  gsl_matrix_float_free(pM2);
  gsl_vector_float_free(pQv);
  gsl_vector_float_free(pV1);
  gsl_vector_float_free(pV2);

  return 0;
}
*/

int getH(EKF_ctx_t *ctx) {
  gsl_matrix_float *H = ctx->wk->H;
  gsl_matrix_float_set_zero(H);
  
  gsl_vector_float_view u_g = gsl_matrix_float_subcolumn(H, 0, 0, 3);
  gsl_vector_float_view u_r = gsl_matrix_float_subcolumn(H, 0, 3, 3);


  gsl_matrix_float_view upperRightH = gsl_matrix_float_submatrix(H, 0, 1, 3, 3);
  gsl_matrix_float_view lowerRightH = gsl_matrix_float_submatrix(H, 3, 1, 3, 3);

  gsl_vector_float *pQv = gsl_vector_float_alloc(3);
  gsl_quat_float_get_imaginary(ctx->q_est, pQv);

  cross_product(ctx->horizonRefG, pQv, &u_g.vector);
  cross_product(ctx->horizonRefMag, pQv, &u_r.vector);

  // upperRightH

  gsl_matrix_float *pM1 = &upperRightH.matrix;
  gsl_matrix_float *pM2 = gsl_matrix_float_calloc(3, 3);
  gsl_vector_float *pV1 = gsl_vector_float_alloc(3);
  gsl_vector_float *pV2 = gsl_vector_float_alloc(3);

  gsl_vector_float_memcpy(pV2, ctx->horizonRefG);
  gsl_vector_float_scale(pV2, gsl_quat_float_get(ctx->q_est, 0));

  gsl_vector_float_memcpy(pV1, &u_g.vector);
  gsl_vector_float_add(pV1, pV2);

  skewSymFromVector(pV1, pM1);

  gsl_blas_sger(1, ctx->horizonRefG, pQv, pM2);
  gsl_matrix_float_scale(pM2, -1.f);
  gsl_matrix_float_add(pM1, pM2);

  gsl_matrix_float_set_identity(pM2);

  float G_Qv_dot = 0;
  gsl_blas_sdot(pQv, ctx->horizonRefG, &G_Qv_dot);
  gsl_matrix_float_scale(pM2, G_Qv_dot);

  gsl_matrix_float_add(pM1, pM2);

  // lowerRightH

  pM1 = &lowerRightH.matrix;
  gsl_matrix_float_set_zero(pM2);

  gsl_vector_float_memcpy(pV2, ctx->horizonRefMag);
  gsl_vector_float_scale(pV2, gsl_quat_float_get(ctx->q_est, 0));

  gsl_vector_float_memcpy(pV1, &u_r.vector);
  gsl_vector_float_add(pV1, pV2);

  skewSymFromVector(pV1, pM1);

  gsl_blas_sger(1, ctx->horizonRefMag, pQv, pM2);
  gsl_matrix_float_scale(pM2, -1.f);
  gsl_matrix_float_add(pM1, pM2);

  gsl_matrix_float_set_identity(pM2);

  float Mag_Qv_dot = 0;
  gsl_blas_sdot(pQv, ctx->horizonRefMag, &Mag_Qv_dot);
  gsl_matrix_float_scale(pM2, Mag_Qv_dot);

  gsl_matrix_float_add(pM1, pM2);
  gsl_matrix_float_scale(H, 2.f);

  gsl_matrix_float_free(pM2);
  gsl_vector_float_free(pQv);
  gsl_vector_float_free(pV1);
  gsl_vector_float_free(pV2);

  return 0;
}



int getR(EKF_ctx_t *ctx) {
  gsl_matrix_float *R = ctx->wk->R;
  gsl_matrix_float_set_zero(R);
  for (int i = 0; i < 3; i++) {
    gsl_matrix_float_set(R, i, i, ACC_STD_DEVIATION);
    gsl_matrix_float_set(R, i + 3, i + 3, MAG_STD_DEVIATION);
  }
  return 0;
}

int getS(EKF_ctx_t *ctx) {
  getR(ctx);
  gsl_matrix_float *S = ctx->wk->S;
  gsl_matrix_float *H = ctx->wk->H;
  gsl_matrix_float_set_zero(ctx->wk->S);

  gsl_matrix_float *bufferMat = gsl_matrix_float_calloc(4, 6);

  gsl_blas_sgemm(CblasNoTrans, CblasTrans, 1, ctx->P_est, H, 0, bufferMat);
  gsl_blas_sgemm(CblasNoTrans, CblasNoTrans, 1, H, bufferMat, 0, S);
  gsl_matrix_float_add(S, ctx->wk->R);

  gsl_matrix_float_free(bufferMat);

  return 0;
}

int getK(EKF_ctx_t *ctx) {
  getH(ctx);
  getS(ctx);

  gsl_matrix_float *K = ctx->wk->K;

  // K = P*H'*inv(S);
  gsl_matrix_float *S = ctx->wk->S;
  gsl_matrix *doubleS = ctx->wk->doubleS;
  gsl_matrix *doubleSLU = gsl_matrix_alloc(doubleS->size1, doubleS->size2);
  for (int i = 0; i < S->size1; i++) {
    for (int j = 0; j < S->size2; j++) {
      gsl_matrix_set(doubleS, i, j, gsl_matrix_float_get(S, i, j));
    }
  }

  gsl_permutation *p = gsl_permutation_alloc(S->size1);
  gsl_matrix_memcpy(doubleSLU, doubleS);
  gsl_matrix *doubleInvS = ctx->wk->doubleInvS;
  gsl_matrix_set_zero(doubleInvS);
  int signum = 0;
  gsl_linalg_LU_decomp(doubleSLU, p, &signum);
  gsl_linalg_LU_invert(doubleSLU, p, doubleInvS);

  gsl_matrix_float *invS = ctx->wk->invS;
  for (int i = 0; i < invS->size1; i++) {
    for (int j = 0; j < invS->size2; j++) {
      gsl_matrix_float_set(invS, i, j, gsl_matrix_get(doubleInvS, i, j));
    }
  }
  gsl_permutation_free(p);

  gsl_blas_sgemm(CblasTrans, CblasNoTrans, 1, ctx->wk->H, invS, 0, ctx->wk->M1_4_6);
  gsl_blas_sgemm(CblasNoTrans, CblasNoTrans, 1, ctx->P_est, ctx->wk->M1_4_6, 0, K);

  gsl_matrix_free(doubleSLU);

  return 0;
}

void ekfNorm(EKF_ctx_t *ctx) { gsl_quat_float_normamilize(ctx->q_current); }

void ekfInitConditions(EKF_ctx_t *ctx, const measures_t *measures) {
  qInitEstimate(ctx, measures);
  PInitEstimate(ctx);
}

void qInitEstimate(EKF_ctx_t *ctx, const measures_t *measures) {
  gsl_vector_float *magMeasure = gsl_vector_float_alloc(3);
  gsl_vector_float *accMeasure = gsl_vector_float_alloc(3);

  gsl_vector_float *xAxis = gsl_vector_float_alloc(3);
  gsl_vector_float *yAxis = gsl_vector_float_alloc(3);
  gsl_vector_float *zAxis = gsl_vector_float_alloc(3);

  gsl_matrix_float *rotMat = gsl_matrix_float_alloc(3, 3);

  float magNorm, accNorm;

  for (uint8_t i = 0; i < 3; i++) {
    gsl_vector_float_set(magMeasure, i, measures->mag[i]);
    gsl_vector_float_set(accMeasure, i, measures->acc[i]);
  }

  gsl_blas_sdot(magMeasure, magMeasure, &magNorm);
  gsl_blas_sdot(accMeasure, accMeasure, &accNorm);

  magNorm = sqrt(magNorm);
  accNorm = sqrt(accNorm);

  gsl_vector_float_scale(magMeasure, 1 / magNorm);
  gsl_vector_float_scale(accMeasure, 1 / accNorm);

  gsl_vector_float_memcpy(zAxis, accMeasure);
  cross_product(zAxis, magMeasure, yAxis);
  cross_product(yAxis, zAxis, xAxis);

  gsl_matrix_float_set_col(rotMat, 0, xAxis);
  gsl_matrix_float_set_col(rotMat, 1, yAxis);
  gsl_matrix_float_set_col(rotMat, 2, zAxis);

  gsl_quat_float *q = gsl_quat_float_alloc();
  gsl_quat_float_fromRotMatrix(rotMat, q);
  gsl_vector_float_memcpy(ctx->q_current, q);
  gsl_quat_float_free(q);

  gsl_vector_float_free(xAxis);
  gsl_vector_float_free(yAxis);
  gsl_vector_float_free(zAxis);
  gsl_vector_float_free(magMeasure);
  gsl_vector_float_free(accMeasure);
}

void PInitEstimate(EKF_ctx_t *ctx) {
  gsl_matrix_float_set_identity(ctx->P_current);
}
