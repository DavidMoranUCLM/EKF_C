add_compile_options(-Wall)

add_executable(TEST_QUATERNION quaternion_test.c)
target_include_directories(TEST_QUATERNION PUBLIC ../include ../unity/src)
target_link_libraries(TEST_QUATERNION unity GLS_QUATERNION)

add_test(quaternion_test TEST_QUATERNION)

add_executable(TEST_JACOBIAN jacobian_test.c)
target_include_directories(TEST_JACOBIAN PUBLIC ../include ../unity/src)
target_link_libraries(TEST_JACOBIAN unity JACOBIAN)

add_test(jacobian_test TEST_JACOBIAN)

add_executable(TEST_EKF EKF_test.c)
target_include_directories(TEST_EKF PUBLIC ../include ../unity/src)
target_link_libraries(TEST_EKF unity GLS_QUATERNION EKF ROTATIONS)

add_test(EKF_test TEST_EKF)