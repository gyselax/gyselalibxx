SET( SPLINES_FILES
  bsplines.cpp
  bsplines_uniform.cpp
  bsplines_non_uniform.cpp
  spline_1d.cpp
  spline_interpolator_1d.cpp
  spline_2d.cpp
  spline_interpolator_2d.cpp
)

add_library( cxx_splines STATIC ${SPLINES_FILES} )
set_property(TARGET cxx_splines PROPERTY Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/cxx_splines_FMOD")
target_include_directories(cxx_splines
    PUBLIC
        "${CMAKE_CURRENT_SOURCE_DIR}"
        "${CMAKE_CURRENT_SOURCE_DIR}/../math_tools"
        "${CMAKE_CURRENT_SOURCE_DIR}/../boundary_conditions"
    INTERFACE
        "${CMAKE_CURRENT_BINARY_DIR}/cxx_splines_FMOD"
    PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries( cxx_splines
  PUBLIC
     std::mdspan
     cxx_spline_matrix
)

add_library(selalib::cxx_splines ALIAS cxx_splines)
