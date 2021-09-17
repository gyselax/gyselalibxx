cmake_minimum_required(VERSION 3.15)
project(splines CXX)

# Our project

## The library itself
add_library(splines STATIC
	src/deprecated/bsplines.cpp
	src/deprecated/bsplines_uniform.cpp
	src/deprecated/bsplines_non_uniform.cpp
	src/deprecated/spline_1d.cpp
	src/deprecated/spline_2d.cpp
	src/deprecated/spline_builder_1d.cpp
	src/deprecated/spline_builder_2d.cpp
	src/matrix.cpp
	src/matrix_dense.cpp
	src/matrix_banded.cpp
	src/matrix_corner_block.cpp
	src/matrix_center_block.cpp
	src/matrix_periodic_banded.cpp
	src/matrix_pds_tridiag.cpp
	src/null_boundary_value.cpp
)
target_compile_features(splines PUBLIC cxx_std_17)
target_include_directories(splines
	PUBLIC
		"$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
)
target_link_libraries(splines
	PUBLIC
		DDC::DDC
	PRIVATE
		${LAPACK_LIBRARIES}
)
add_library(sll::splines ALIAS splines)

## if tests are enabled, build the tests in `tests/`
if("${BUILD_TESTING}")
    add_subdirectory(tests/)
endif()
