#include "EKF.h"

#include "EKF_const.h"
#include "gsl_quaternion_float.h"
#include "gsl/gsl_math.h"
#include "rotations.h"
#include "unity.h"

#define FLOAT_ERROR 1e-6f
#define ESTIMATE_ERROR 1e-1f
#define CORRECTION_ERROR 1e-2f
#define TIME_STEP 1e-1f

/**
 * Private declarations
 *
 **/

EKF_ctx_t EKF_ctx;
EKF_work_ctx_t EKF_wk_ctx;

/**
 * Private function declarations
 *
 **/
void setUp(void);
void tearDown(void);

void suiteSetUp(void);

void resetTest(void);
void verifyTest(void);

void testInvFloat(void);
void testGet_h(void);
void testGetH(void);
void testEKFInit(void);
void testStep(void);

/**
 * Main
 *
 **/

int main() {
  UNITY_BEGIN();
  // RUN_TEST(testEKFInit);
  RUN_TEST(testInvFloat);
  RUN_TEST(testGet_h);
  RUN_TEST(testGetH);
  RUN_TEST(testStep);
  return UNITY_END();
}

/**
 * Private function definitions
 *
 **/

void setUp(void) {
  float acc[3] = {0, 0, 9.8};
  float mag[3] = {cosf(LATITUDE_RAD), 0, -sinf(LATITUDE_RAD)};
  float velAng[3] = {0, 0, 0};
  measures_t measure;
  measure.acc[0] = acc[0];
  measure.acc[1] = acc[1];
  measure.acc[2] = acc[2];
  measure.mag[0] = mag[0];
  measure.mag[1] = mag[1];
  measure.mag[2] = mag[2];
  measure.velAng[0] = velAng[0];
  measure.velAng[1] = velAng[1];
  measure.velAng[2] = velAng[2];

  ekfInit(&EKF_ctx, &EKF_wk_ctx, &measure);
}
void tearDown(void) { ekfDeinit(&EKF_ctx); }

void suiteSetUp(void) {}
// int suiteTearDown(int num_failures) {}

void resetTest(void) {}
void verifyTest(void) {}



void testInvFloat(void){
  gsl_matrix_float *S = gsl_matrix_float_calloc(6,6);
  gsl_matrix_float *invS = gsl_matrix_float_alloc(6,6);
  gsl_matrix_float_set_identity(S);
  invertMatrixFloat(&EKF_ctx, S, invS);

  TEST_ASSERT_FLOAT_ARRAY_WITHIN(FLOAT_ERROR,S->data, invS->data, 6*6);

  gsl_matrix_float_set_zero(invS);

  gsl_matrix_float_set_identity(S);
  gsl_matrix_float_scale(S, 4);

  invertMatrixFloat(&EKF_ctx, S, invS);

  gsl_matrix_float_set_identity(S);
  gsl_matrix_float_scale(S, 0.25);

  TEST_ASSERT_FLOAT_ARRAY_WITHIN(FLOAT_ERROR, S->data, invS->data, 6*6);
}

void testEKFInit(void) {
  float q0[4] = {1, 0, 0, 0};
  float acc[3] = {0, 0, 9.8};
  float mag[3] = {1, 0, 0};
  float velAng[3] = {0, 0, 0};
  measures_t measure;
  measure.acc[0] = acc[0];
  measure.acc[1] = acc[1];
  measure.acc[2] = acc[2];
  measure.mag[0] = mag[0];
  measure.mag[1] = mag[1];
  measure.mag[2] = mag[2];
  measure.velAng[0] = velAng[0];
  measure.velAng[1] = velAng[1];
  measure.velAng[2] = velAng[2];

  ekfInit(&EKF_ctx, &EKF_wk_ctx, &measure);
  TEST_ASSERT_FLOAT_ARRAY_WITHIN(FLOAT_ERROR, acc, EKF_ctx.acc->data, 3);
  TEST_ASSERT_FLOAT_ARRAY_WITHIN(FLOAT_ERROR, mag, EKF_ctx.mag->data, 3);
  TEST_ASSERT_FLOAT_ARRAY_WITHIN(FLOAT_ERROR, velAng, EKF_ctx.velAng->data, 3);

  TEST_ASSERT_EQUAL_FLOAT_ARRAY(q0, EKF_ctx.q_current->data, 4);

  gsl_matrix_float *P0 = gsl_matrix_float_calloc(4, 4);
  gsl_matrix_float_set_identity(P0);
  TEST_ASSERT_EQUAL_FLOAT_ARRAY(P0->data, EKF_ctx.P_current->data, 4 * 4);

  ekfDeinit(&EKF_ctx);
}

void testGet_h(void) {
  float rotAngRad = 3.14 / 2;
  gsl_vector_float *Axis = gsl_vector_float_calloc(3);
  gsl_vector_float_set(Axis, 0, 0);
  gsl_vector_float_set(Axis, 1, 0);
  gsl_vector_float_set(Axis, 2, 1);
  gsl_quat_float_fromAxis(Axis, rotAngRad, EKF_ctx.q_est);
  get_h(&EKF_ctx);
  TEST_ASSERT_FLOAT_ARRAY_WITHIN(FLOAT_ERROR, EKF_ctx.horizonRefG->data,
                                 EKF_ctx.wk->h->data, 3);
  TEST_ASSERT_FLOAT_ARRAY_WITHIN(FLOAT_ERROR, EKF_ctx.horizonRefMag->data,
                                 EKF_ctx.wk->h->data + 3, 3);
}

// 0	0	2	0
// 0	-2	0	0
// 0	0	0	0
// 0	0	-1.28	0
// 0	1.28	0	1.53
// 0	0	-1.53	0

void testGetH(void) {
  float rotAngRad = 20.0*3.1415/180.0;
  gsl_vector_float *Axis = gsl_vector_float_calloc(3);
  gsl_vector_float_set(Axis, 0, 0);
  gsl_vector_float_set(Axis, 1, 0);
  gsl_vector_float_set(Axis, 2, 1);
  gsl_quat_float_fromAxis(Axis, rotAngRad, EKF_ctx.q_est);
  getH(&EKF_ctx);
}

void testStep(void) {
  //float q0[4] = {1, 0, 0, 0};

  gsl_vector_float *pAcc = gsl_vector_float_calloc(3);
  gsl_vector_float *pMag = gsl_vector_float_calloc(3);
  gsl_vector_float *pVelAng = gsl_vector_float_calloc(3);

  gsl_vector_float_set(pAcc, 0, 0);
  gsl_vector_float_set(pAcc, 1, 0);
  gsl_vector_float_set(pAcc, 2, 9.8);

  gsl_vector_float_set(pMag, 0, cosf(EKF_ctx.latitude));
  gsl_vector_float_set(pMag, 1, 0);
  gsl_vector_float_set(pMag, 2, -sinf(EKF_ctx.latitude));

  gsl_vector_float_set(pVelAng, 0, 0);
  gsl_vector_float_set(pVelAng, 1, 0);
  gsl_vector_float_set(pVelAng, 2, 0.1);

  measures_t measure;

  gsl_quat_float *pQw = gsl_quat_float_alloc();
  gsl_vector_float *pV = gsl_vector_float_calloc(3);
  gsl_vector_float_set(pV, 2, 1);
  gsl_quat_float_fromAxis(pV, -0.1 * TIME_STEP, pQw);

  rotation_t rotation;
  createRotationFromQuat(pQw, &rotation);
  for (int i = 0; i < 100000; i++) {
    

    measure.acc[0] = gsl_vector_float_get(pAcc, 0);
    measure.acc[1] = gsl_vector_float_get(pAcc, 1);
    measure.acc[2] = gsl_vector_float_get(pAcc, 2);
    measure.mag[0] = gsl_vector_float_get(pMag, 0);
    measure.mag[1] = gsl_vector_float_get(pMag, 1);
    measure.mag[2] = gsl_vector_float_get(pMag, 2);
    measure.velAng[0] = gsl_vector_float_get(pVelAng, 0);
    measure.velAng[1] = gsl_vector_float_get(pVelAng, 1);
    measure.velAng[2] = gsl_vector_float_get(pVelAng, 2);

    ekfStep(&EKF_ctx, &measure, TIME_STEP*i);


    gsl_quat_float_fromAxis(pV, -0.1 * TIME_STEP*i, pQw);
    createRotationFromQuat(pQw, &rotation);

    gsl_vector_float_set(pAcc, 0, 0);
    gsl_vector_float_set(pAcc, 1, 0);
    gsl_vector_float_set(pAcc, 2, 9.8);

    gsl_vector_float_set(pMag, 0, cosf(EKF_ctx.latitude));
    gsl_vector_float_set(pMag, 1, 0);
    gsl_vector_float_set(pMag, 2, -sinf(EKF_ctx.latitude));

    rotateVector(pAcc, &rotation);
    rotateVector(pMag, &rotation);

    // TEST_ASSERT_FLOAT_ARRAY_WITHIN(ESTIMATE_ERROR, q0, EKF_ctx.q_est->data,
    // 4); TEST_ASSERT_FLOAT_ARRAY_WITHIN(ESTIMATE_ERROR, q0,
    // EKF_ctx.q_current->data,
    //                                4);
  }
}