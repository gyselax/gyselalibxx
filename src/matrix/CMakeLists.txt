add_library(matrix STATIC
	matrix.cpp
	matrix_dense.cpp
	matrix_banded.cpp
	matrix_corner_block.cpp
	matrix_center_block.cpp
	matrix_periodic_banded.cpp
	matrix_pds_tridiag.cpp
)

target_compile_features(matrix
	PUBLIC
		cxx_std_17
)

target_include_directories(matrix
	PUBLIC
		"${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries(matrix
	PUBLIC
		std::mdspan
		sll::math_tools
	PRIVATE
		${LAPACK_LIBRARIES}
)

add_library(sll::matrix ALIAS matrix)
