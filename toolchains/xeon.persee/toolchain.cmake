
include(${CMAKE_CURRENT_LIST_DIR}/../common_toolchains/importable_defaults.cmake)

# CMake options
set(CMAKE_BUILD_TYPE Release)

# Compiler options
set(CMAKE_CXX_FLAGS "-Wall -Wno-sign-compare -Wno-unused-but-set-variable")

# Kokkos options
set(Kokkos_ENABLE_OPENMP ON CACHE BOOL "Activate OpenMP parallelisation" FORCE)

# Activate/deactivate parts of the code
set(ACTIVATE_RESTART_TESTS OFF)
