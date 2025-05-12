#include <inttypes.h>

#define ACC_STD_DEVIATION (0.5*0.5)
#define GYRO_STD_DEVIATION (0.15*0.15)
#define MAG_STD_DEVIATION (0.5*0.5)
#define GYRO_BIAS_NOISE (0.005*0.005)

#define STATE_SIZE 7
#define CORRECTION_SIZE 3

#define ESTIMATE_ORDER (uint8_t)1
#define DIFERENTIAL_STEP 1e-3F
#define P_ESTIMATE_SCALE 1.F

#define MAG_CORRECTION_PERIOD_S 1.f


#define ACC_SCALE 1.F
#define MAG_SCALE 1.F

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
