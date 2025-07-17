
include(${CMAKE_CURRENT_LIST_DIR}/../common_toolchains/importable_defaults.cmake)

# CMake options
set(CMAKE_BUILD_TYPE Release)

# Compiler options
set(CMAKE_CXX_FLAGS "-Wall -Wno-sign-compare -Wno-unused-but-set-variable")

# Kokkos options
set(Kokkos_ENABLE_OPENMP ON CACHE BOOL "Activate OpenMP parallelisation" FORCE)

# Kokkos Kernels options
set(KokkosKernels_ENABLE_ALL_COMPONENTS OFF CACHE BOOL "")
set(KokkosKernels_ENABLE_COMPONENT_BATCHED ON CACHE BOOL "")
set(KokkosKernels_ENABLE_COMPONENT_BLAS ON CACHE BOOL "")
set(KokkosKernels_ADD_DEFAULT_ETI OFF CACHE BOOL "")

# Activate/deactivate parts of the code
set(ACTIVATE_RESTART_TESTS OFF)
