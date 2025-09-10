
include(${CMAKE_CURRENT_LIST_DIR}/../../common_toolchains/importable_defaults.cmake)

# CMake options
set(CMAKE_BUILD_TYPE Release)

# Compiler options
set(CMAKE_CXX_COMPILER nvcc_wrapper)
set(CMAKE_CXX_EXTENSIONS OFF) # Avoid a Kokkos warning that will force if to OFF anyway when compiling with nvcc
# The compile option ipa-sra triggers a segfault with nvcc. @tpadioleau reported it to Nvidia. We then disable it.
# Using CUDA 12.9 triggers a spurious warning when compiling for Nvidia Volta architectures. We then disable it.
set(CMAKE_CXX_FLAGS_INIT "-fno-ipa-sra -Wno-deprecated-gpu-targets")
set(CMAKE_CXX_FLAGS_INIT "${CMAKE_CXX_FLAGS_INIT} -Wall -Wno-sign-compare --Werror cross-execution-space-call -Xcudafe --diag_suppress=unsigned_compare_with_zero -Xcudafe --diag_suppress=integer_sign_change")

# Gyselalibxx options
set(GYSELALIBXX_DEFAULT_CXX_FLAGS "" CACHE STRING "Default flags for C++ specific to Gyselalib++" FORCE)

# Activate/deactivate parts of the code
if (DEFINED ENV{DDC_BUILD_TESTING})
    set(DDC_BUILD_TESTS $ENV{DDC_BUILD_TESTING})
endif()
set(ACTIVATE_RESTART_TESTS OFF)
