
include(${CMAKE_CURRENT_LIST_DIR}/../cicd_default_toolchain.cmake)

# CMake options
set(CMAKE_BUILD_TYPE FakeDebug)

# Compiler options
set(CMAKE_CXX_FLAGS "-Wall -Werror=parentheses -Werror=unused-local-typedefs -Werror=unused-variable -Wno-sign-compare")

# Kokkos options
set(Kokkos_ENABLE_SERIAL ON CACHE BOOL "Allow serial code to run" FORCE)

# Activate/deactivate parts of the code
set(BUILD_BENCHMARKS ON)
set(ACTIVATE_RESTART_TESTS ON)
