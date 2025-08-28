set(CMAKE_CXX_COMPILER hipcc)
set(CMAKE_C_COMPILER amdclang)
set(CMAKE_Fortran_COMPILER amdflang)
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
set(CMAKE_CXX_FLAGS_INIT "${CMAKE_CXX_FLAGS_INIT} -isystem $ENV{SPACK_USER_PREFIX}/linux-rhel8-zen3/gcc-13.2.1.mi250.raw/openblas-0.3.28-ipbz/include")
set(CMAKE_EXE_LINKER_FLAGS_INIT "$ENV{PE_MPICH_GTL_DIR_amd_gfx90a} $ENV{PE_MPICH_GTL_LIBS_amd_gfx90a}")
