
include(${CMAKE_CURRENT_LIST_DIR}/../cicd_default_toolchain.cmake)

# CMake options
set(CMAKE_BUILD_TYPE Release)

# Compiler options
set(CMAKE_CXX_FLAGS "-Wall -Wno-sign-compare -Wno-unused-but-set-variable")

# Activate/deactivate parts of the code
if (DEFINED ENV{DDC_BUILD_TESTING})
    set(DDC_BUILD_TESTS $ENV{DDC_BUILD_TESTING})
endif()
set(ACTIVATE_RESTART_TESTS OFF)
