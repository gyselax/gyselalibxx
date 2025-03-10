set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/installation")

set(CMAKE_CXX_COMPILER hipcc)
set(CMAKE_C_COMPILER amdclang)
set(CMAKE_Fortran_COMPILER amdflang)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
# NOTE: We are not supposed to define CMAKE_BUILD_TYPE here.
set(CMAKE_BUILD_TYPE Release) # Debug, Release, RelWithDebInfo and MinSizeRel

# Kokkos options:

set(Kokkos_ENABLE_HIP ON CACHE BOOL "Defined if the HIP execution space is enabled.")
set(Kokkos_ENABLE_HIP_MULTIPLE_KERNEL_INSTANTIATIONS ON CACHE BOOL "If defined, multiple kernel versions are instantiated potentially improving run time.")
set(Kokkos_ARCH_VEGA90A ON CACHE BOOL "Enable support for AMD GPU MI200 series (GFX90A).")
# FIXME: According to Thomas this is needed.In practice I did not help the
# build. Note that it can be detrimental to performance.
set(Kokkos_ENABLE_HIP_RELOCATABLE_DEVICE_CODE ON CACHE BOOL "")

set(Kokkos_ENABLE_OPENMP OFF CACHE BOOL "")
set(Kokkos_ENABLE_SERIAL ON CACHE BOOL "")
set(Kokkos_ARCH_ZEN3 ON CACHE BOOL "Optimize for AMD Zen3 architecture (HOST).")

# Kokkos kernels options:

set(KokkosKernels_ENABLE_ALL_COMPONENTS OFF CACHE BOOL "")
set(KokkosKernels_ENABLE_COMPONENT_BATCHED ON CACHE BOOL "")
set(KokkosKernels_ENABLE_COMPONENT_BLAS ON CACHE BOOL "")

# Koliop options:

set(koliop_ENABLE_LTO OFF CACHE BOOL "")

# The rest is optional:

# We should add that too, but there is too much warnings ! -Wsuggest-override -Wctor-dtor-privacy -Wdouble-promotion -Wcast-qual -Wredundant-decls -Wswitch-default -Wold-style-cast -Wswitch-enum -Wundef
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic -Wcast-align -Wformat=2 -Winit-self -Woverloaded-virtual -Wsign-promo -Wstrict-aliasing -Wdisabled-optimization -Wtautological-compare -Wpacked -Wunreachable-code -Wno-sign-compare -Wno-unused-parameter -Wno-unused-but-set-variable")

# FIXME: blas are not properly detected.
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem $ENV{CRAY_LIBSCI_PREFIX}/include" CACHE STRING "")

# NOTE: People often export this as environment variable instead.
set(BLAS_LIBRARIES "$ENV{CRAY_LIBSCI_PREFIX}/lib/libsci_cray.so;" CACHE STRING "")
# set(LAPACK_LIBRARIES "$ENV{CRAY_LIBSCI_PREFIX}/lib/libsci_cray.so;" CACHE STRING "")
# set(LAPACKE_LIBRARIES "$ENV{CRAY_LIBSCI_PREFIX}/lib/libsci_cray.so;" CACHE STRING "")

# FIXME: SLL compiler crash
set(SLL_BUILD_TESTING OFF CACHE BOOL "")
