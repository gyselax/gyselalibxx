add_library(cxx_spline_matrix STATIC
  matrix.cpp
  matrix_dense.cpp
  matrix_banded.cpp
  matrix_corner_block.cpp
  matrix_center_block.cpp
  matrix_periodic_banded.cpp
  matrix_pds_tridiag.cpp
 )
target_include_directories(cxx_spline_matrix
    PUBLIC
        "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries(cxx_spline_matrix
  PUBLIC
    std::mdspan
    sll::math_tools
  PRIVATE
    ${LAPACK_LIBRARIES}
)

add_library(sll::cxx_spline_matrix ALIAS cxx_spline_matrix)
