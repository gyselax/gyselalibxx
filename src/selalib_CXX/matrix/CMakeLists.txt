SET( MATRIX_CXX_FILES
  matrix.cpp
  matrix_dense.cpp
  matrix_banded.cpp
  matrix_corner_block.cpp
  matrix_center_block.cpp
  matrix_periodic_banded.cpp
  matrix_pds_tridiag.cpp
  )

add_library( cxx_spline_matrix OBJECT ${MATRIX_CXX_FILES} )
target_include_directories(cxx_spline_matrix
    PUBLIC
        "${CMAKE_CURRENT_SOURCE_DIR}"
    PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}/../math_tools"
)

target_link_libraries( cxx_spline_matrix
  PUBLIC
    std::mdspan
  PRIVATE
    ${LAPACK_LIBRARIES}
)

add_library(selalib::cxx_spline_matrix ALIAS cxx_spline_matrix)
