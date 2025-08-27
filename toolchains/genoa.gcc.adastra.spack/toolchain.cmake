set(CMAKE_CXX_COMPILER g++)
set(CMAKE_C_COMPILER gcc)
set(CMAKE_Fortran_COMPILER gfortran)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
# NOTE: We are not supposed to define CMAKE_BUILD_TYPE here.
set(CMAKE_BUILD_TYPE Release) # Debug, Release, RelWithDebInfo and MinSizeRel

# Koliop options:

set(koliop_ENABLE_LTO OFF CACHE BOOL "")

# The rest is optional:

# We should add that too, but there is too much warnings ! -Wsuggest-override -Wctor-dtor-privacy -Wdouble-promotion -Wcast-qual -Wredundant-decls -Wswitch-default -Wold-style-cast -Wswitch-enum -Wundef
set(CMAKE_CXX_FLAGS_INIT "${CMAKE_CXX_FLAGS_INIT} -Wall -Wextra -Wpedantic -Wcast-align -Wformat=2 -Winit-self -Woverloaded-virtual -Wsign-promo -Wstrict-aliasing -Wdisabled-optimization -Wtautological-compare -Wpacked -Wunreachable-code -Wno-sign-compare -Wno-unused-parameter -Wno-unused-but-set-variable")

# FIXME: blas are not properly detected.
set(CMAKE_CXX_FLAGS_INIT "${CMAKE_CXX_FLAGS_INIT} -isystem $ENV{CRAY_LIBSCI_PREFIX}/include")

# NOTE: People often export this as environment variable instead.
set(BLAS_LIBRARIES "$ENV{CRAY_LIBSCI_PREFIX}/lib/libsci_gnu_mp.so;" CACHE STRING "")
# set(LAPACK_LIBRARIES "$ENV{CRAY_LIBSCI_PREFIX}/lib/libsci_gnu_mp.so;" CACHE STRING "")
# set(LAPACKE_LIBRARIES "$ENV{CRAY_LIBSCI_PREFIX}/lib/libsci_gnu_mp.so;" CACHE STRING "")
