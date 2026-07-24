
include(${CMAKE_CURRENT_LIST_DIR}/importable_defaults.cmake)

# CMake options
set(CMAKE_BUILD_TYPE Debug)

# Compiler options
set(CMAKE_CXX_FLAGS_INIT "-O1 -Wall -Wno-sign-compare")
set(CMAKE_COMPILE_WARNING_AS_ERROR ON)

# Activate/deactivate parts of the code
set(GYSELALIBXX_ACTIVATE_RESTART_TESTS ON)
