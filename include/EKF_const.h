#include <inttypes.h>

#define ACC_STD_DEVIATION (0.4*0.4)
#define MAG_STD_DEVIATION (0.3*0.3)
#define GYRO_STD_DEVIATION (0.2*0.2)

#define P_SIZE 4

#define ESTIMATE_ORDER (uint8_t)1
#define DIFERENTIAL_STEP 1e-3F

#define LATITUDE_DEG 40.F
#define LATITUDE_RAD ((LATITUDE_DEG*3.1415/180.f)) 
