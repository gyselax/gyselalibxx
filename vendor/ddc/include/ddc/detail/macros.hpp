#pragma once

#include <Kokkos_Core.hpp>

#if defined(__NVCC__)

#define DDC_MAKE_PRAGMA(X) _Pragma(#X)

#if defined __NVCC_DIAG_PRAGMA_SUPPORT__
#define DDC_NV_DIAGNOSTIC_PUSH DDC_MAKE_PRAGMA(nv_diagnostic push)
#define DDC_NV_DIAGNOSTIC_POP DDC_MAKE_PRAGMA(nv_diagnostic pop)
#define DDC_NV_DIAG_SUPPRESS(X) DDC_MAKE_PRAGMA(nv_diag_suppress X)
#else
#define DDC_NV_DIAGNOSTIC_PUSH
#define DDC_NV_DIAGNOSTIC_POP
#define DDC_NV_DIAG_SUPPRESS(X) DDC_MAKE_PRAGMA(diag_suppress X)
#endif

#define DDC_IF_NVCC_THEN_PUSH(X) DDC_NV_DIAGNOSTIC_PUSH
#define DDC_IF_NVCC_THEN_SUPPRESS(X) DDC_NV_DIAG_SUPPRESS(X)
#define DDC_IF_NVCC_THEN_PUSH_AND_SUPPRESS(X)                                                      \
    DDC_NV_DIAGNOSTIC_PUSH                                                                         \
    DDC_NV_DIAG_SUPPRESS(X)
#define DDC_IF_NVCC_THEN_POP DDC_NV_DIAGNOSTIC_POP

#else

#define DDC_IF_NVCC_THEN_PUSH(X)
#define DDC_IF_NVCC_THEN_SUPPRESS(X)
#define DDC_IF_NVCC_THEN_PUSH_AND_SUPPRESS(X)
#define DDC_IF_NVCC_THEN_POP

#endif

#if defined(__HIP_DEVICE_COMPILE__)
#define DDC_CURRENT_KOKKOS_SPACE Kokkos::HIPSpace
#elif defined(__CUDA_ARCH__)
#define DDC_CURRENT_KOKKOS_SPACE Kokkos::CudaSpace
#else
#define DDC_CURRENT_KOKKOS_SPACE Kokkos::HostSpace
#endif
