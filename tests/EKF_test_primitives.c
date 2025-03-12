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
#include "string.h"
#include "unity.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// Input file names
#define M3_FILE   "/media/david/data/Users/deivi/Documents/Asignaturas/TFG/test_primitives/M3.bin"
#define M4_FILE   "/media/david/data/Users/deivi/Documents/Asignaturas/TFG/test_primitives/M4.bin"
#define Q4_FILE   "/media/david/data/Users/deivi/Documents/Asignaturas/TFG/test_primitives/Q4.bin"
#define M6_FILE   "/media/david/data/Users/deivi/Documents/Asignaturas/TFG/test_primitives/M6.bin"
#define M3_3_FILE "/media/david/data/Users/deivi/Documents/Asignaturas/TFG/test_primitives/M3_3.bin"
#define M4_4_FILE "/media/david/data/Users/deivi/Documents/Asignaturas/TFG/test_primitives/M4_4.bin"
#define M4_6_FILE "/media/david/data/Users/deivi/Documents/Asignaturas/TFG/test_primitives/M4_6.bin"
#define M6_4_FILE "/media/david/data/Users/deivi/Documents/Asignaturas/TFG/test_primitives/M6_4.bin"
#define M6_6_FILE "/media/david/data/Users/deivi/Documents/Asignaturas/TFG/test_primitives/M6_6.bin"

// Base output file names
#define Q_EST_OUT_FILE_BASE     "/media/david/data/Users/deivi/Documents/Asignaturas/TFG/test_primitives/q_est_out_C.bin"
#define P_EST_OUT_FILE_BASE     "/media/david/data/Users/deivi/Documents/Asignaturas/TFG/test_primitives/P_est_out_C.bin"
#define P_CORRECT_OUT_FILE_BASE "/media/david/data/Users/deivi/Documents/Asignaturas/TFG/test_primitives/P_correct_out_C.bin"
#define Q_CORRECT_OUT_FILE_BASE "/media/david/data/Users/deivi/Documents/Asignaturas/TFG/test_primitives/q_correct_out_C.bin"
#define K_OUT_FILE_BASE         "/media/david/data/Users/deivi/Documents/Asignaturas/TFG/test_primitives/K_out_C.bin"
#define H_OUT_FILE_BASE         "/media/david/data/Users/deivi/Documents/Asignaturas/TFG/test_primitives/H_out_C.bin"
#define h_OUT_FILE_BASE         "/media/david/data/Users/deivi/Documents/Asignaturas/TFG/test_primitives/h_out_C.bin"
#define R_OUT_FILE_BASE         "/media/david/data/Users/deivi/Documents/Asignaturas/TFG/test_primitives/R_out_C.bin"
#define S_OUT_FILE_BASE         "/media/david/data/Users/deivi/Documents/Asignaturas/TFG/test_primitives/S_out_C.bin"
#define INV_S_OUT_FILE_BASE "/media/david/data/Users/deivi/Documents/Asignaturas/TFG/test_primitives/invS_out_C.bin"

#define N_SAMPLES 100 // Define the number of samples (n)

// Helper function to read a gsl_vector_float from a binary file, with an index
gsl_vector_float *read_vector(const char *filename, int size, int index) {
  FILE *f = fopen(filename, "rb");
  if (!f) {
    printf("Failed to open file: %s\n", filename);
    return NULL;
  }
  gsl_vector_float *v = gsl_vector_float_alloc(size);
  fseek(f, index * size * sizeof(float), SEEK_SET); // Seek to the correct position
  size_t read_count = fread(v->data, sizeof(float), size, f);
  fclose(f);
  if (read_count != size) {
    printf("Failed to read vector from %s (read %zu elements, expected %d)\n",
           filename, read_count, size);
    gsl_vector_float_free(v);
    return NULL;
  }
  return v;
}

// Helper function to read a gsl_matrix_float from a binary file, with an index
gsl_matrix_float *read_matrix(const char *filename, int size1, int size2, int index) {
  FILE *f = fopen(filename, "rb");
  if (!f) {
    printf("Failed to open file: %s\n", filename);
    return NULL;
  }
  gsl_matrix_float *m = gsl_matrix_float_alloc(size1, size2);
  fseek(f, index * size1 * size2 * sizeof(float), SEEK_SET); // Seek to the correct position
  size_t read_count = fread(m->data, sizeof(float), size1 * size2, f);
  fclose(f);
  if (read_count != size1 * size2) {
    printf("Failed to read matrix from %s (read %zu elements, expected %d)\n",
           filename, read_count, size1 * size2);
    gsl_matrix_float_free(m);
    return NULL;
  }
  return m;
}

// Helper function to write a gsl_vector_float to a binary file, with a filename
void write_vector_to_file(const char *filename, size_t index, const gsl_vector_float *v) {
  FILE *f = fopen(filename, "ab");
  if (f) {
    fseek(f, index * v->size * sizeof(float), SEEK_SET);
    fwrite(v->data, sizeof(float), v->size, f);
    fclose(f);
  } else {
    printf("Failed to open file for writing: %s\n", filename);
  }
}

// Helper function to write a gsl_matrix_float to a binary file, with a filename
void write_matrix_to_file(const char *filename, size_t index, const gsl_matrix_float *m) {
  FILE *f = fopen(filename, "ab");
  if (f) {
    fseek(f, index * m->size1 * m->size2 * sizeof(float), SEEK_SET);
    fwrite(m->data, sizeof(float), m->size1 * m->size2, f);
    fclose(f);
  } else {
    printf("Failed to open file for writing: %s\n", filename);
  }
}

void setUp(void) {
  // set stuff up here
}

void tearDown(void) {
  // clean stuff up here
}

void test_qEstPrimitive(void) {
  for (int i = 0; i < N_SAMPLES; i++) {
    gsl_vector_float *velAng = read_vector(M3_FILE, 3, i);
    gsl_quat_float *qPrev = read_vector(Q4_FILE, 4, i);
    gsl_quat_float *qEst = gsl_quat_float_alloc();
    gsl_quat_float *tmpQuat = gsl_quat_float_alloc();
    gsl_matrix_float *qVelAngMat = gsl_matrix_float_alloc(4, 4);
    gsl_matrix_float *q1Mat = gsl_matrix_float_alloc(4, 4);
    float deltaT = 0.01f;

    qEstPrimitive(velAng, deltaT, qPrev, qEst, tmpQuat, qVelAngMat, q1Mat);
    
    write_vector_to_file(Q_EST_OUT_FILE_BASE, i, (gsl_vector_float *)qEst);

    gsl_vector_float_free(velAng);
    gsl_quat_float_free(qPrev);
    gsl_quat_float_free(qEst);
    gsl_quat_float_free(tmpQuat);
    gsl_matrix_float_free(qVelAngMat);
    gsl_matrix_float_free(q1Mat);
  }
}

void test_hPrimitive(void) {
  for (int i = 0; i < N_SAMPLES; i++) {
    gsl_matrix_float *PPrev = read_matrix(M4_4_FILE, 4, 4, i);
    gsl_matrix_float *F = read_matrix(M4_4_FILE, 4, 4, i);
    gsl_matrix_float *Q = read_matrix(M4_4_FILE, 4, 4, i);
    gsl_matrix_float *PEst = gsl_matrix_float_alloc(4, 4);
    gsl_matrix_float *tmp4x4 = gsl_matrix_float_alloc(4, 4);

    PEstPrimitive(PPrev, F, Q, PEst, tmp4x4);

    write_matrix_to_file(P_EST_OUT_FILE_BASE, i, PEst);

    gsl_matrix_float_free(PPrev);
    gsl_matrix_float_free(F);
    gsl_matrix_float_free(Q);
    gsl_matrix_float_free(PEst);
    gsl_matrix_float_free(tmp4x4);
  }
}


void test_PEstPrimitive(void) {
  for (int i = 0; i < N_SAMPLES; i++) {
    gsl_quat_float *q = read_vector(Q4_FILE, 4, i);
    gsl_vector_float *r = read_vector(M3_FILE, 3, i);
    gsl_vector_float *g = read_vector(M3_FILE, 3, i);
    gsl_vector_float *h = gsl_vector_float_calloc(6);

    get_hPrimitive(q, g, r, h);

    write_vector_to_file(h_OUT_FILE_BASE, i, h);

    gsl_vector_float_free(q);
    gsl_vector_float_free(r);
    gsl_vector_float_free(g);
    gsl_vector_float_free(h);
  }
}



void test_PCorrectPrimitive(void) {
  for (int i = 0; i < N_SAMPLES; i++) {
    gsl_matrix_float *P = read_matrix(M4_4_FILE, 4, 4, i);
    gsl_matrix_float *K = read_matrix(M4_6_FILE, 4, 6, i);
    gsl_matrix_float *H = read_matrix(M6_4_FILE, 6, 4, i);
    gsl_matrix_float *R = read_matrix(M6_6_FILE, 6, 6, i);
    gsl_matrix_float *PCorrect = gsl_matrix_float_alloc(4, 4);
    gsl_matrix_float *tmp4x4 = gsl_matrix_float_alloc(4, 4);
    gsl_matrix_float *tmp6x4 = gsl_matrix_float_alloc(6, 4);
    gsl_matrix_float *I4 = gsl_matrix_float_alloc(4, 4);

    PCorrectPrimitive(P, K, H, R, PCorrect, tmp4x4, tmp6x4, I4);

    write_matrix_to_file(P_CORRECT_OUT_FILE_BASE, i, PCorrect);

    gsl_matrix_float_free(P);
    gsl_matrix_float_free(K);
    gsl_matrix_float_free(H);
    gsl_matrix_float_free(R);
    gsl_matrix_float_free(PCorrect);
    gsl_matrix_float_free(tmp4x4);
    gsl_matrix_float_free(tmp6x4);
    gsl_matrix_float_free(I4);
  }
}

void test_qCorrectPrimitive(void) {
  for (int i = 0; i < N_SAMPLES; i++) {
    gsl_quat_float *q = read_vector(Q4_FILE, 4, i);
    gsl_matrix_float *K = read_matrix(M4_6_FILE, 4, 6, i);
    gsl_vector_float *h = read_vector(M6_FILE, 6, i);
    gsl_quat_float *qCorrect = gsl_quat_float_alloc();
    gsl_vector_float *tmp = gsl_vector_float_alloc(4);

    qCorrectPrimitive(q, K, h, qCorrect, tmp);

    write_vector_to_file(Q_CORRECT_OUT_FILE_BASE, i, (gsl_vector_float *)qCorrect);

    gsl_quat_float_free(q);
    gsl_matrix_float_free(K);
    gsl_vector_float_free(h);
    gsl_quat_float_free(qCorrect);
    gsl_vector_float_free(tmp);
  }
}

void test_getKPrimitive(void) {
  for (int i = 0; i < N_SAMPLES; i++) {
    gsl_matrix_float *P = read_matrix(M4_4_FILE, 4, 4, i);
    gsl_matrix_float *H = read_matrix(M6_4_FILE, 6, 4, i);
    gsl_matrix_float *invS = read_matrix(M6_6_FILE, 6, 6, i);
    gsl_matrix_float *K = gsl_matrix_float_alloc(4, 6);
    gsl_matrix_float *tmp4x6 = gsl_matrix_float_alloc(4, 6);

    getKPrimitive(P, H, invS, K, tmp4x6);

    write_matrix_to_file(K_OUT_FILE_BASE, i, K);

    gsl_matrix_float_free(P);
    gsl_matrix_float_free(H);
    gsl_matrix_float_free(invS);
    gsl_matrix_float_free(K);
    gsl_matrix_float_free(tmp4x6);
  }
}

void test_getHPrimitive(void) {
  for (int i = 0; i < N_SAMPLES; i++) {
    gsl_quat_float *q = read_vector(Q4_FILE, 4, i);
    gsl_vector_float *acc = read_vector(M3_FILE, 3, i);
    gsl_vector_float *mag = read_vector(M3_FILE, 3, i);
    gsl_vector_float *velAng = read_vector(M3_FILE, 3, i);
    gsl_matrix_float *H = gsl_matrix_float_alloc(6, 4);
    gsl_vector_float *pQv = gsl_vector_float_alloc(3);
    gsl_matrix_float *pM2 = gsl_matrix_float_alloc(3, 3);
    gsl_vector_float *pV1 = gsl_vector_float_alloc(3);
    gsl_vector_float *pV2 = gsl_vector_float_alloc(3);

    getHPrimitive(q, acc, mag, H, pQv, pM2, pV1, pV2);

    write_matrix_to_file(H_OUT_FILE_BASE, i, H);

    gsl_quat_float_free(q);
    gsl_vector_float_free(acc);
    gsl_vector_float_free(mag);
    gsl_vector_float_free(velAng);
    gsl_matrix_float_free(H);
    gsl_vector_float_free(pQv);
    gsl_matrix_float_free(pM2);
    gsl_vector_float_free(pV1);
    gsl_vector_float_free(pV2);
  }
}

void test_getRPrimitive(void) {
  for (int i = 0; i < N_SAMPLES; i++) {
    gsl_matrix_float *R = gsl_matrix_float_alloc(6, 6);

    getRPrimitive(R);

    write_matrix_to_file(R_OUT_FILE_BASE, i, R);

    gsl_matrix_float_free(R);
  }
}

void test_getSPrimitive(void) {
  for (int i = 0; i < N_SAMPLES; i++) {
    gsl_matrix_float *H = read_matrix(M6_4_FILE, 6, 4, i);
    gsl_matrix_float *P = read_matrix(M4_4_FILE, 4, 4, i);
    gsl_matrix_float *R = read_matrix(M6_6_FILE, 6, 6, i);
    gsl_matrix_float *S = gsl_matrix_float_alloc(6, 6);
    gsl_matrix_float *tmp4x6 = gsl_matrix_float_alloc(4, 6);

    getSPrimitive(H, P, R, S, tmp4x6);

    write_matrix_to_file(S_OUT_FILE_BASE, i, S);

    gsl_matrix_float_free(H);
    gsl_matrix_float_free(P);
    gsl_matrix_float_free(R);
    gsl_matrix_float_free(S);
    gsl_matrix_float_free(tmp4x6);
  }
}

void test_invertMatrixFloat(void) {
    for (int i = 0; i < N_SAMPLES; i++) {
        gsl_matrix_float *S = read_matrix(M6_6_FILE, 6, 6, i);
        gsl_matrix_float *invS = gsl_matrix_float_alloc(6, 6);
        gsl_matrix *tmpMatrix_d = gsl_matrix_alloc(6, 6);
        gsl_matrix *doubleInvS = gsl_matrix_alloc(6, 6);
        gsl_vector *tmpTau_d = gsl_vector_alloc(6);
        gsl_vector *tmpB_d = gsl_vector_alloc(6);
        gsl_vector *tmpX_d = gsl_vector_alloc(6);
        
        int n = S->size1;
        // Convert float matrix S into double-precision matrix using preallocated buffer
        for (int k = 0; k < n; k++) {
            for (int j = 0; j < n; j++) {
                gsl_matrix_set(tmpMatrix_d, k, j, gsl_matrix_float_get(S, k, j));
            }
        }
        invertMatrixFloatPrimitive(tmpMatrix_d, doubleInvS, tmpTau_d, tmpB_d, tmpX_d);
        for (int k = 0; k < n; k++) {
            for (int j = 0; j < n; j++) {
                gsl_matrix_float_set(invS, k, j, (float)gsl_matrix_get(doubleInvS, k, j));
            }
        }
        write_matrix_to_file(INV_S_OUT_FILE_BASE, 0, invS);

        gsl_matrix_float_free(S);
        gsl_matrix_float_free(invS);
        gsl_matrix_free(tmpMatrix_d);
        gsl_matrix_free(doubleInvS);
        gsl_vector_free(tmpTau_d);
        gsl_vector_free(tmpB_d);
        gsl_vector_free(tmpX_d);
    }
}

int main(void) {
  UNITY_BEGIN();
  RUN_TEST(test_qEstPrimitive);
  RUN_TEST(test_PEstPrimitive);
  RUN_TEST(test_PCorrectPrimitive);
  RUN_TEST(test_qCorrectPrimitive);
  RUN_TEST(test_hPrimitive);
  RUN_TEST(test_getKPrimitive);
  RUN_TEST(test_getHPrimitive);
  RUN_TEST(test_getRPrimitive);
  RUN_TEST(test_getSPrimitive);
  RUN_TEST(test_invertMatrixFloat);
  return UNITY_END();
}