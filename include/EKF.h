#ifndef EKF_H
#define EKF_H

#include "gsl/gsl_matrix.h"
#include "gsl/gsl_matrix_float.h"
#include "gsl/gsl_vector_float.h"
#include "gsl_quaternion_float.h"

typedef struct measures_s {
  float acc[3];
  float velAng[3];
  float mag[3];
} measures_t;

typedef struct EKF_work_ctx_s {
  // Estimation
  gsl_matrix_float *F;
  gsl_matrix_float *W;
  gsl_matrix_float *Q;
  // Correction
  gsl_matrix_float *K;
  gsl_matrix_float *H;
  gsl_matrix_float *S;
  gsl_matrix_float *invS;
  gsl_matrix *doubleS;
  gsl_matrix *doubleInvS;
  gsl_matrix_float *R;

  gsl_vector_float *z;
  gsl_vector_float *h;

  // Buffers temporales (existing)
  const gsl_matrix_float *I3;
  const gsl_matrix_float *I4;
  gsl_matrix_float *M1_4_3;
  gsl_matrix_float *M2_4_4;
  gsl_vector_float *v1;
  gsl_vector_float *v2;

  // NEW: Preallocated temporary buffers to avoid alloc/free on every call  
  gsl_quat_float  *tmpQuat;         // 4-element quaternion buffer  
  gsl_matrix_float *tmp4x4;          // 4x4 matrix buffer  
  gsl_matrix_float *tmp3x3;          // 3x3 matrix buffer  
  gsl_matrix_float *tmpStdDevMat;    // 3x3 temporary for standard deviation  
  gsl_matrix_float *tmpBufferMat;    // 3x4 temporary buffer  
  gsl_vector_float *tmp3vec;         // 3-element vector temporary  

  // NEW: Additional preallocated double-precision buffers for QR inversion
  gsl_matrix      *inv_tmpMatrix_d; // double matrix of size 6x6
  gsl_vector      *inv_tmpTau_d;    // double vector of size 6
  gsl_vector      *inv_tmpB_d;      // double vector of size 6
  gsl_vector      *inv_tmpX_d;      // double vector of size 6

  // NEW: Temporary matrix for the product R * Trans(K)
  gsl_matrix_float *tmpRTransK;

} EKF_work_ctx_t;

typedef struct EKF_ctx_s {
  EKF_work_ctx_t *wk;
  gsl_quat_float *q_current;
  gsl_quat_float *q_prev;
  gsl_quat_float *q_est;

  gsl_matrix_float *P_prev;
  gsl_matrix_float *P_current;
  gsl_matrix_float *P_est;

  gsl_vector_float *acc;
  gsl_vector_float *velAng;
  gsl_vector_float *mag;

  gsl_vector_float *magStdDev;

  gsl_vector_float *horizonRefG;

  float lastMagCorrectionTime_s;
  float magCorrectionPeriod_s;
  float currentTime;
  float prevTime;
  float lastMagCorection;

  int (*getSemaphore)(void* sem);
  int (*releaseSemaphore)(void* sem);

} EKF_ctx_t;

/**
 * @brief Initializes a EKF_ctx_t pointer in heap with a reference measurement
 *
 * @param ctx
 * @param measures
 */
void ekfInit(EKF_ctx_t *ctx, const measures_t *measures);

/**
 * @brief Deinitializes a EKF_ctx_t pointer, freeing memory;
 *
 * @param ctx
 */
void ekfDeinit(EKF_ctx_t *ctx);

/**
 * @brief Updates the filter with new measuments of a given time
 *
 * @param ctx
 * @param measures
 * @param currentTime
 */
void ekfStep(EKF_ctx_t *ctx, const measures_t *measures,
             const float currentTime);


#ifdef TEST_ENABLE

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
void invertMatrixFloat(EKF_ctx_t *ctx, const gsl_matrix_float *S,
                       gsl_matrix_float *invS);
int getK(EKF_ctx_t *ctx);
void ekfNorm(EKF_ctx_t *ctx);
void ekfInitConditions(EKF_ctx_t *ctx, const measures_t *measures);
void qInitEstimate(EKF_ctx_t *ctx, const measures_t *measures);
void PInitEstimate(EKF_ctx_t *ctx);

void qEstPrimitive(const gsl_vector_float *velAng, float deltaT,
                   const gsl_quat_float *qPrev, gsl_quat_float *qEst,
                   gsl_quat_float *tmpQuat, gsl_matrix_float *qVelAngMat,
                   gsl_matrix_float *q1Mat);
void PEstPrimitive(const gsl_matrix_float *PPrev, const gsl_matrix_float *F,
                   const gsl_matrix_float *Q, gsl_matrix_float *PEst,
                   gsl_matrix_float *tmp4x4);
void PCorrectPrimitive(const gsl_matrix_float *P, const gsl_matrix_float *K,
                       const gsl_matrix_float *H, const gsl_matrix_float *R,
                       gsl_matrix_float *PCorrect, gsl_matrix_float *tmp4x4,
                       gsl_matrix_float *tmp3x4, gsl_matrix_float *I4);
void qCorrectPrimitive(const gsl_quat_float *q, const gsl_matrix_float *K,
                       const gsl_vector_float *h, gsl_quat_float *qCorrect,
                       gsl_vector_float *tmp);
void getKPrimitive(const gsl_matrix_float *P, const gsl_matrix_float *H,
                   const gsl_matrix_float *invS, gsl_matrix_float *K,
                   gsl_matrix_float *tmp4x3);
void getHPrimitive(const gsl_quat_float *q, const gsl_vector_float *acc,
                   gsl_matrix_float *H, gsl_vector_float *pQv,
                   gsl_matrix_float *pM2, gsl_vector_float *pV1, gsl_vector_float *pV2);
void getRPrimitive(gsl_matrix_float *R);
void getSPrimitive(const gsl_matrix_float *H, const gsl_matrix_float *P,
                   const gsl_matrix_float *R, gsl_matrix_float *S,
                   gsl_matrix_float *tmp4x3);
void invertMatrixFloatPrimitive(const gsl_matrix *S, gsl_matrix *invS,
                                gsl_vector *tau, gsl_vector *b, gsl_vector *x);
void get_hPrimitive(const gsl_quat_float *q_est,
                    const gsl_vector_float *horizonRefG, gsl_vector_float *h);
#endif

#endif