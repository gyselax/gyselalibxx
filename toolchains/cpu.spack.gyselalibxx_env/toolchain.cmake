
include(${CMAKE_CURRENT_LIST_DIR}/../common_toolchains/importable_defaults.cmake)

# Specify Python version to ensure spack packages are used
set(Python3_EXECUTABLE $ENV{PYTHON_EXECUTABLE})

# CMake options
set(CMAKE_BUILD_TYPE Release)

# Kokkos Kernels options
set(KokkosKernels_ENABLE_ALL_COMPONENTS OFF CACHE BOOL "")
set(KokkosKernels_ENABLE_COMPONENT_BATCHED ON CACHE BOOL "")
set(KokkosKernels_ENABLE_COMPONENT_BLAS ON CACHE BOOL "")
set(KokkosKernels_ADD_DEFAULT_ETI OFF CACHE BOOL "")

# Compiler options
set(CMAKE_CXX_FLAGS_INIT "-Wall -Wno-sign-compare")

# Activate/deactivate parts of the code
set(BUILD_BENCHMARKS ON)
set(ACTIVATE_RESTART_TESTS ON)
