#pragma once

#include <Kokkos_Core.hpp>

#if defined(__NVCC__)
#define DDC_INTERNAL_FIX_NVCC_IF_CONSTEXPR

#define DDC_MAKE_PRAGMA(X) _Pragma(#X)

#if defined __NVCC_DIAG_PRAGMA_SUPPORT__
#define DDC_NV_DIAG_SUPPRESS(X) DDC_MAKE_PRAGMA(nv_diag_suppress X)
#define DDC_NV_DIAG_DEFAULT(X) DDC_MAKE_PRAGMA(nv_diag_default X)
#else
#define DDC_NV_DIAG_SUPPRESS(X) DDC_MAKE_PRAGMA(diag_suppress X)
#define DDC_NV_DIAG_DEFAULT(X) DDC_MAKE_PRAGMA(diag_default X)
#endif

#endif

#define DDC_LAMBDA KOKKOS_LAMBDA

#define DDC_INLINE_FUNCTION KOKKOS_INLINE_FUNCTION

#define DDC_FORCEINLINE_FUNCTION KOKKOS_FORCEINLINE_FUNCTION
