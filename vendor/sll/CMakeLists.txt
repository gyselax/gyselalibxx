cmake_minimum_required(VERSION 3.15)
project(SLL CXX)

# Our project
option(SLL_BUILD_TESTING "Build sll tests" ON)

## The library itself
add_library(SLL STATIC
    src/matrix.cpp
    src/matrix_dense.cpp
    src/matrix_banded.cpp
    src/matrix_corner_block.cpp
    src/matrix_center_block.cpp
    src/matrix_periodic_banded.cpp
    src/matrix_pds_tridiag.cpp
    src/gauss_legendre_integration.cpp
)
target_compile_features(SLL PUBLIC cxx_std_17)
target_include_directories(SLL
    PUBLIC
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
)
target_link_libraries(SLL
    PUBLIC
        DDC::core
        DDC::splines
				Kokkos::kokkos
)
if(CMAKE_VERSION GREATER_EQUAL 3.18)
    target_link_libraries(SLL PRIVATE LAPACK::LAPACK)
else()
    target_link_libraries(SLL PRIVATE ${LAPACK_LIBRARIES})
endif()
add_library(sll::SLL ALIAS SLL)

## if tests are enabled, build the tests in `tests/`
if("${GYSELALIBXX_BUILD_TESTING}" AND "${SLL_BUILD_TESTING}")
    add_subdirectory(tests/)
endif()
