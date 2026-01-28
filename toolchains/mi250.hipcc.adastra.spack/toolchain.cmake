# CMake options
# NOTE: We are not supposed to define CMAKE_BUILD_TYPE here.
set(CMAKE_BUILD_TYPE Release) # Debug, Release, RelWithDebInfo and MinSizeRel

# Compiler options
set(CMAKE_CXX_COMPILER hipcc)
set(CMAKE_C_COMPILER amdclang)
set(CMAKE_Fortran_COMPILER amdflang)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
# We should add that too, but there is too much warnings ! -Wsuggest-override -Wctor-dtor-privacy -Wdouble-promotion -Wcast-qual -Wredundant-decls -Wswitch-default -Wold-style-cast -Wswitch-enum -Wundef
set(CMAKE_CXX_FLAGS_INIT "-Wall -Wextra -Wpedantic -Wcast-align -Wformat=2 -Winit-self -Woverloaded-virtual -Wsign-promo -Wstrict-aliasing -Wdisabled-optimization -Wtautological-compare -Wpacked -Wunreachable-code -Wno-sign-compare -Wno-unused-parameter -Wno-unused-but-set-variable")
set(CMAKE_CXX_FLAGS_INIT "${CMAKE_CXX_FLAGS_INIT} -isystem $ENV{GYSELALIBXX_OPENBLAS_ROOT}/include")
set(CMAKE_EXE_LINKER_FLAGS_INIT "-L/opt/cray/pe/mpich/8.1.30/gtl/lib -lmpi_gtl_hsa")
