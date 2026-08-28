
include(${CMAKE_CURRENT_LIST_DIR}/importable_defaults.cmake)

# CMake options
set(CMAKE_BUILD_TYPE Release)

# Compiler options
set(CMAKE_CXX_FLAGS_INIT "-Wall -Wextra -Wno-sign-compare -Wno-unused-but-set-parameter -Wno-unused-parameter")

# Activate/deactivate parts of the code
set(GYSELALIBXX_ACTIVATE_RESTART_TESTS OFF)
