
include(${CMAKE_CURRENT_LIST_DIR}/../cicd_default_toolchain.cmake)

# CMake options
set(CMAKE_BUILD_TYPE Debug)

# Compiler options
set(CMAKE_CXX_COMPILER nvcc_wrapper)
set(CMAKE_CXX_FLAGS "-g -Wall -Werror -Wno-sign-compare -Xcudafe --diag_suppress=unsigned_compare_with_zero -Xcudafe --diag_suppress=integer_sign_change -Wno-unused-but-set-variable")

# Activate/deactivate parts of the code
if (DEFINED ENV{DDC_BUILD_TESTING})
    set(DDC_BUILD_TESTS $ENV{DDC_BUILD_TESTING})
endif()
set(ACTIVATE_RESTART_TESTS OFF)
