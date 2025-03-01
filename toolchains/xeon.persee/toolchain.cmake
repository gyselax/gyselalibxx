
include(${CMAKE_CURRENT_LIST_DIR}/../cicd_default_toolchain.cmake)

# CMake options
set(CMAKE_BUILD_TYPE Release)

# Compiler options
set(CMAKE_CXX_FLAGS "-Wall -Wno-sign-compare -Wno-unused-but-set-variable")

# Activate/deactivate parts of the code
set(ACTIVATE_RESTART_TESTS OFF)
