
include(${CMAKE_CURRENT_LIST_DIR}/importable_defaults.cmake)

# CMake options
set(CMAKE_BUILD_TYPE Release)

# Compiler options
set(CMAKE_CXX_FLAGS "-Wall -Werror=unused-local-typedefs -Werror=unused-variable -Wno-sign-compare")

# Kokkos options
set(Kokkos_ENABLE_SERIAL ON CACHE BOOL "Allow serial code to run" FORCE)

# Activate/deactivate parts of the code
set(BUILD_BENCHMARKS ON)
set(ACTIVATE_RESTART_TESTS ON)
