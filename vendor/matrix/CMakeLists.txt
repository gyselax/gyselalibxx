cmake_minimum_required(VERSION 3.15)

add_library(matrix STATIC
	src/matrix.cpp
	src/matrix_dense.cpp
	src/matrix_banded.cpp
	src/matrix_corner_block.cpp
	src/matrix_center_block.cpp
	src/matrix_periodic_banded.cpp
	src/matrix_pds_tridiag.cpp
)

target_compile_features(matrix PUBLIC cxx_std_17)

target_include_directories(matrix
	PUBLIC
		"$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
)

target_link_libraries(matrix
	PUBLIC
		mdspan_wrapper::mdsw
	PRIVATE
		${LAPACK_LIBRARIES}
)

add_library(sll::matrix ALIAS matrix)

