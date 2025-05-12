/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: utils_rotmatYawCorrectionJac.h
 *
 * MATLAB Coder version            : 24.2
 * C/C++ source code generated on  : 22-Apr-2025 21:43:47
 */

#ifndef UTILS_YAWCORRECTIONJAC_H
#define UTILS_YAWCORRECTIONJAC_H

/* Include Files */
// #include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
extern void rotmatYawCorrectionJac(const float in1[7], const float in2[3],
                                   float qJ[7][3]);
extern void quatYawCorrectionJac(const float in1[7], const float in2[3],
                                 float qJ[7][3]);

extern void quatYawCorrectionJacS(const float in1[7], const float in2[3],
                                 float qJS[7][7]);
#ifdef __cplusplus
}
#endif

#endif
/*
 * File trailer for utils_rotmatYawCorrectionJac.h
 *
 * [EOF]
 */
