#include <inttypes.h>

#define ACC_STD_DEVIATION (0.8*0.8)
#define MAG_STD_DEVIATION (1.4*1.4)
#define GYRO_STD_DEVIATION (0.5*0.5)

#define P_SIZE 4

#define ESTIMATE_ORDER (uint8_t)1
#define DIFERENTIAL_STEP 1e-3F

#define LATITUDE_DEG 5.F
#define LATITUDE_RAD ((LATITUDE_DEG*3.1415/180.f)) 
