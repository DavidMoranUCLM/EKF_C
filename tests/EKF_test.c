#include "EKF.h"

#include "EKF_const.h"
#include "gsl/gsl_math.h"
#include "gsl/gsl_matrix_float.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_rng.h"
#include "gsl/gsl_vector_float.h"
#include "gsl_quaternion_float.h"
#include "rotations.h"
#include "stdio.h"
#include "unity.h"

#define FLOAT_ERROR 1e-6f
#define ESTIMATE_ERROR 1e-1f
#define CORRECTION_ERROR 1e-2f
#define TIME_STEP 1e-2f
#define MAX_LOGS_NUMBER 50

#define NUM_ITERATIONS 50000

/**
 * Private declarations
 *
 **/

EKF_ctx_t EKF_ctx;
EKF_work_ctx_t EKF_wk_ctx;

typedef enum matrix_vector_e {
  GSL_VECTOR,
  GSL_MATRIX,
} matrix_vector_t;

typedef struct csvLog_s {
  uint64_t index;
  matrix_vector_t type;
  void *pMV;
  char filename[50];
  void *buffer;
} csvLog_t;

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

void logMatrixCSV_init(csvLog_t *csv, void *pM, const char filename[],
                       matrix_vector_t type);
void logMatrixCSV_update(csvLog_t *csv);
void logMatrixCSV_updateAll();
void logMatrixCSV_deinit(csvLog_t *csv);
void logMatrixCSV_deinitAll();

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

void testInvFloat(void) {
  gsl_matrix_float *S = gsl_matrix_float_calloc(6, 6);
  gsl_matrix_float *invS = gsl_matrix_float_alloc(6, 6);
  gsl_matrix_float_set_identity(S);
  invertMatrixFloat(&EKF_ctx, S, invS);

  TEST_ASSERT_FLOAT_ARRAY_WITHIN(FLOAT_ERROR, S->data, invS->data, 6 * 6);

  gsl_matrix_float_set_zero(invS);

  gsl_matrix_float_set_identity(S);
  gsl_matrix_float_scale(S, 4);

  invertMatrixFloat(&EKF_ctx, S, invS);

  gsl_matrix_float_set_identity(S);
  gsl_matrix_float_scale(S, 0.25);

  TEST_ASSERT_FLOAT_ARRAY_WITHIN(FLOAT_ERROR, S->data, invS->data, 6 * 6);
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
  float rotAngRad = 20.0 * 3.1415 / 180.0;
  gsl_vector_float *Axis = gsl_vector_float_calloc(3);
  gsl_vector_float_set(Axis, 0, 0);
  gsl_vector_float_set(Axis, 1, 0);
  gsl_vector_float_set(Axis, 2, 1);
  gsl_quat_float_fromAxis(Axis, rotAngRad, EKF_ctx.q_est);
  getH(&EKF_ctx);
}

void testStep(void) {
  // float q0[4] = {1, 0, 0, 0};

  gsl_vector_float *pAcc = gsl_vector_float_calloc(3);
  gsl_vector_float *pMag = gsl_vector_float_calloc(3);
  gsl_vector_float *pVelAng = gsl_vector_float_calloc(3);

  gsl_vector_float_set(pAcc, 0, 0);
  gsl_vector_float_set(pAcc, 1, 0);
  gsl_vector_float_set(pAcc, 2, 9.8);

  gsl_vector_float_set(pMag, 0, cosf(EKF_ctx.latitude));
  gsl_vector_float_set(pMag, 1, 0);
  gsl_vector_float_set(pMag, 2, -sinf(EKF_ctx.latitude));

  gsl_vector_float_set(pVelAng, 0, 0.1);
  gsl_vector_float_set(pVelAng, 1, 0);
  gsl_vector_float_set(pVelAng, 2, 0);

  float VelAngNorm;
  gsl_blas_sdot(pVelAng, pVelAng, &VelAngNorm);
  VelAngNorm = sqrtf(VelAngNorm);

  measures_t measure;

  gsl_quat_float *pQw = gsl_quat_float_alloc();
  gsl_quat_float_fromAxis(pVelAng, 0, pQw);

  rotation_t rotation;
  createRotationFromQuat(pQw, &rotation);

  csvLog_t magLog, accLog, quatLog, expectedQuatLog, vLog, HLog, PLog, FLog,
      SLog, PestLog, qEstLog, WLog, QLog, RLog, invSLog, KLog;
  logMatrixCSV_init(&magLog, pMag, "magLog.txt", GSL_VECTOR);
  logMatrixCSV_init(&accLog, pAcc, "accLog.txt", GSL_VECTOR);
  logMatrixCSV_init(&quatLog, EKF_ctx.q_current, "quatLog.txt", GSL_VECTOR);
  logMatrixCSV_init(&qEstLog, EKF_ctx.q_est, "qEstLog.txt", GSL_VECTOR);
  logMatrixCSV_init(&expectedQuatLog, pQw, "quatExpectedLog.txt", GSL_VECTOR);
  logMatrixCSV_init(&vLog, EKF_ctx.wk->z, "vLog.txt", GSL_VECTOR);
  logMatrixCSV_init(&HLog, EKF_ctx.wk->H, "HLog.txt", GSL_MATRIX);
  logMatrixCSV_init(&PLog, EKF_ctx.P_current, "PLog.txt", GSL_MATRIX);
  logMatrixCSV_init(&PestLog, EKF_ctx.P_est, "PestLog.txt", GSL_MATRIX);
  logMatrixCSV_init(&FLog, EKF_ctx.wk->F, "FLog.txt", GSL_MATRIX);
  logMatrixCSV_init(&SLog, EKF_ctx.wk->S, "SLog.txt", GSL_MATRIX);
  logMatrixCSV_init(&invSLog, EKF_ctx.wk->invS, "invSLog.txt", GSL_MATRIX);
  logMatrixCSV_init(&WLog, EKF_ctx.wk->W, "WLog.txt", GSL_MATRIX);
  logMatrixCSV_init(&QLog, EKF_ctx.wk->Q, "QLog.txt", GSL_MATRIX);
  logMatrixCSV_init(&RLog, EKF_ctx.wk->R, "RLog.txt", GSL_MATRIX);
  logMatrixCSV_init(&KLog, EKF_ctx.wk->K, "KLog.txt", GSL_MATRIX);

  gsl_rng *randHandle = gsl_rng_alloc(gsl_rng_default);

  for (int i = 0; i < NUM_ITERATIONS; i++) {
    measure.acc[0] = gsl_vector_float_get(pAcc, 0);
    measure.acc[1] = gsl_vector_float_get(pAcc, 1);
    measure.acc[2] = gsl_vector_float_get(pAcc, 2);
    measure.mag[0] = gsl_vector_float_get(pMag, 0);
    measure.mag[1] = gsl_vector_float_get(pMag, 1);
    measure.mag[2] = gsl_vector_float_get(pMag, 2);
    measure.velAng[0] = gsl_vector_float_get(pVelAng, 0);
    measure.velAng[1] = gsl_vector_float_get(pVelAng, 1);
    measure.velAng[2] = gsl_vector_float_get(pVelAng, 2);

    ekfStep(&EKF_ctx, &measure, TIME_STEP * i);

    gsl_quat_float_fromAxis(pVelAng, VelAngNorm * TIME_STEP * i, pQw);
    gsl_quat_float_conjugate(
        pQw);  // Giro de sistema de coordenadas, no de vector
    createRotationFromQuat(pQw, &rotation);

    gsl_vector_float_set(pAcc, 0, 0 + (float)gsl_ran_gaussian(randHandle, 0));
    gsl_vector_float_set(pAcc, 1, 0 + (float)gsl_ran_gaussian(randHandle, 0));
    gsl_vector_float_set(pAcc, 2,
                         9.8 + (float)gsl_ran_gaussian(randHandle, 0));

    gsl_vector_float_set(pMag, 0, cosf(EKF_ctx.latitude));
    gsl_vector_float_set(pMag, 1, 0);
    gsl_vector_float_set(pMag, 2, -sinf(EKF_ctx.latitude));

    rotateVector(pMag, &rotation);
    rotateVector(pAcc, &rotation);

    gsl_quat_float_conjugate(pQw);
    logMatrixCSV_updateAll();
  }
  logMatrixCSV_deinitAll();
}

csvLog_t *allCSVLogs[MAX_LOGS_NUMBER] = {NULL};
void logMatrixCSV_init(csvLog_t *csv, void *pMV, const char filename[],
                       matrix_vector_t type) {
  csv->index = 0;
  snprintf(csv->filename, sizeof(csv->filename), "%s", filename);
  csv->type = type;
  csv->pMV = pMV;

  switch (type) {
    case GSL_VECTOR:
      gsl_vector_float *pV = pMV;
      csv->buffer = malloc(pV->size * sizeof(float) * NUM_ITERATIONS);
      break;
    case GSL_MATRIX:
      gsl_matrix_float *pM = pMV;
      csv->buffer =
          malloc(pM->size1 * pM->size2 * sizeof(float) * NUM_ITERATIONS);
    default:
      break;
  }

  for (int i = 0; i < MAX_LOGS_NUMBER; i++) {
    if (allCSVLogs[i] == NULL) {
      allCSVLogs[i] = csv;
      break;
    }
  }
}

void logMatrixCSV_update(csvLog_t *csv) {
  float *buffer = csv->buffer;

  switch (csv->type) {
    case GSL_VECTOR:
      gsl_vector_float *pV = csv->pMV;

      for (int i = 0; i <= pV->size - 1; i++) {
        buffer[(csv->index) * pV->size + i] = gsl_vector_float_get(pV, i);
      }

      break;
    case GSL_MATRIX:
      gsl_matrix_float *pM = csv->pMV;

      for (int i = 0; i <= pM->size1 - 1; i++) {
        for (int j = 0; j <= pM->size2 - 1; j++) {
          buffer[(csv->index) * pM->size1 * pM->size2 + i+j] =
              gsl_matrix_float_get(pM, i, j);
        }
      }

    default:
      break;
  }

  csv->index++;
}

void logMatrixCSV_updateAll() {
  for (int i = 0; i < MAX_LOGS_NUMBER; i++) {
    if (allCSVLogs[i] != NULL) {
      logMatrixCSV_update(allCSVLogs[i]);
    }
  }
}

void logMatrixCSV_deinit(csvLog_t *csv) {
  FILE *file = fopen(csv->filename, "w");
  float *buffer = csv->buffer;

  switch (csv->type) {
    case GSL_VECTOR:
      gsl_vector_float *pV = csv->pMV;
      for (int i = 0; i <= pV->size - 2; i++) {
        fprintf(file, "a%d,", i);
      }
      fprintf(file, "a%lu\n", pV->size - 1);

      for (int row = 0; row < NUM_ITERATIONS; row++) {
        for (int i = 0; i <= pV->size - 2; i++) {
          fprintf(file, "%f,", buffer[row * pV->size + i]);
        }
        fprintf(file, "%f\n", buffer[row * pV->size - 1]);
      }

      break;
    case GSL_MATRIX:
      gsl_matrix_float *pM = csv->pMV;

      for (int i = 0; i <= pM->size1 - 1; i++) {
        for (int j = 0; j <= pM->size2 - 1; j++) {
          fprintf(file, "a%d-%d,", i, j);
        }
      }
      fprintf(file, "a%lu-%lu\n", pM->size1 - 1, pM->size2 - 1);

      for (int row = 0; row < NUM_ITERATIONS; row++) {
        for (int i = 0; i <= pM->size1 - 1; i++) {
          for (int j = 0; j <= pM->size2 - 1; j++) {
            fprintf(file, "%f,", buffer[row * pM->size1 * pM->size2 + i + j]);
          }
        }
        fprintf(file, "%f\n",
                buffer[row * pM->size1 * pM->size2 - 1]);
      }

    default:
      break;
  }

  free(buffer);
}

void logMatrixCSV_deinitAll() {
  for (int i = 0; i < MAX_LOGS_NUMBER; i++) {
    if (allCSVLogs[i] != NULL) {
      logMatrixCSV_deinit(allCSVLogs[i]);
      allCSVLogs[i] = NULL;
    }
  }
}