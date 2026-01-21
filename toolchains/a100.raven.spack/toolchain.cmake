
# CMake options
# NOTE: We are not supposed to define CMAKE_BUILD_TYPE here.
set(CMAKE_BUILD_TYPE Release) # Debug, Release, RelWithDebInfo and MinSizeRel

# Compiler options
set(CMAKE_CXX_COMPILER nvcc_wrapper)
set(CMAKE_CXX_EXTENSIONS OFF) # Avoid a Kokkos warning that will force if to OFF anyway when compiling with nvcc
# The compile option ipa-sra triggers a segfault with nvcc. @tpadioleau reported it to Nvidia. We then disable it.
set(CMAKE_CXX_FLAGS_INIT "-fno-ipa-sra")
set(CMAKE_CXX_FLAGS_INIT "${CMAKE_CXX_FLAGS_INIT} -Wall -Wno-sign-compare --Werror cross-execution-space-call -Xcudafe --diag_suppress=unsigned_compare_with_zero -Xcudafe --diag_suppress=integer_sign_change")
set(CMAKE_CXX_FLAGS_INIT "${CMAKE_CXX_FLAGS_INIT} -isystem $ENV{GYSELALIBXX_OPENBLAS_ROOT}/include")

# Gyselalibxx options
set(GYSELALIBXX_DEFAULT_CXX_FLAGS "" CACHE STRING "Default flags for C++ specific to Gyselalib++" FORCE)

# Koliop options
set(koliop_ENABLE_LTO OFF CACHE BOOL "")
