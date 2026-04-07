
include(${CMAKE_CURRENT_LIST_DIR}/importable_defaults.cmake)

# CMake options
set(CMAKE_BUILD_TYPE Debug)

# Compiler options
set(CMAKE_CXX_FLAGS_INIT "-Og --coverage -fprofile-update=atomic")

# Gyselalibxx options
set(GYSELALIBXX_DEFAULT_CXX_FLAGS "" CACHE STRING "Default flags for C++ specific to Gyselalib++" FORCE)

# Activate/deactivate parts of the code
set(GYSELALIBXX_ACTIVATE_RESTART_TESTS ON)
