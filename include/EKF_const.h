#include <inttypes.h>
#include "gsl/gsl_const.h"

#define ACC_STD_DEVIATION (5*5)
#define MAG_STD_DEVIATION (0.05*0.05)
#define GYRO_STD_DEVIATION (0.15*0.15)

#define P_SIZE 4

#define ESTIMATE_ORDER (uint8_t)1
#define DIFERENTIAL_STEP 1e-3F
#define P_ESTIMATE_SCALE 1.F

#define LATITUDE_DEG 25.F
#define INCLINATION_DEG 0.F
#define LATITUDE_RAD ((LATITUDE_DEG*3.1415f/180.f)) 
#define INCLINATION_RAD ((INCLINATION_DEG*3.1415f/180.f))

#define MAG_SCALE 1.F
#define ACC_SCALE GSL_CONST_MKS_GRAV_ACCEL

#define MAG_CORRECTION_PERIOD_S 5.F

//0: Standard
//1: Joseph
#define P_CORRECT_METHOD 1

/*
#include <inttypes.h>

#define ACC_STD_DEVIATION (0.5*0.5)
#define MAG_STD_DEVIATION (0.8*0.8)
#define GYRO_STD_DEVIATION (0.3*0.3)

#define P_SIZE 4

#define ESTIMATE_ORDER (uint8_t)1
#define DIFERENTIAL_STEP 1e-3F

#define LATITUDE_DEG 50.F
#define LATITUDE_RAD ((LATITUDE_DEG*3.1415/180.f))  
 */