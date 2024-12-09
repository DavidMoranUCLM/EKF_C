cmake_minimum_required(VERSION 3.5.0)
project(INS_C VERSION 0.1.0 LANGUAGES C)

set(GDBINIT_FILE "${CMAKE_SOURCE_DIR}/.gdbinit")

add_compile_definitions(TEST_ENABLE)

include(CTest)
enable_testing()

add_subdirectory(gsl)
link_libraries(gsl gslcblas m)

add_subdirectory(src)


if(NOT DEFINED NOT_TESTING)

    add_subdirectory(tests)
    add_subdirectory(unity)
    
endif()



# Install the headers
install(
    DIRECTORY include/
    DESTINATION include
)

# Install the library
install(
    TARGETS EKF GLS_QUATERNION MATH_UTILS ROTATIONS
    ARCHIVE DESTINATION lib
)
