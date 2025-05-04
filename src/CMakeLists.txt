add_compile_options(-Wall)


add_library(EKF STATIC EKF.c)
add_library(GLS_QUATERNION STATIC gsl_quaternion_float.c)
add_library(JACOBIAN STATIC jacobian.c)
add_library(MATH_UTILS math_utils.c)
add_library(ROTATIONS rotations.c)
add_library(MAG_CORRECTION mag_correction.c)
add_library(YAW_CORRECTION yawCorrectionJac.c)

target_include_directories(EKF PUBLIC ../include ../gsl/vector ../gsl/matrix)
target_include_directories(GLS_QUATERNION PUBLIC ../include ../gsl/vector ../gsl/matrix)
target_include_directories(JACOBIAN PUBLIC ../include ../gsl/vector)
target_include_directories(MATH_UTILS PUBLIC ../include ../gsl/vector)
target_include_directories(ROTATIONS PUBLIC ../include ../gsl/vector ../gsl/matrix)
target_include_directories(MAG_CORRECTION PUBLIC ../include ../gsl/vector ../gsl/matrix)
target_include_directories(YAW_CORRECTION PUBLIC ../include)

target_link_libraries(YAW_CORRECTION GLS_QUATERNION)
target_link_libraries(EKF GLS_QUATERNION MATH_UTILS ROTATIONS MAG_CORRECTION YAW_CORRECTION)