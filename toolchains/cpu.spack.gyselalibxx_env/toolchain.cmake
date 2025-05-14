
include(${CMAKE_CURRENT_LIST_DIR}/../cicd_default_toolchain.cmake)

# Specify Python version to ensure spack packages are used
set(Python3_EXECUTABLE <SPACK_PYTHON_PATH>)

# CMake options
set(CMAKE_BUILD_TYPE Release)

# Compiler options
set(CMAKE_CXX_FLAGS "-Wall -Wno-sign-compare")

# Activate/deactivate parts of the code
set(BUILD_BENCHMARKS ON)
set(ACTIVATE_RESTART_TESTS ON)
