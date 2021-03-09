add_library( cxx_splines STATIC
  bsplines.cpp
  bsplines_uniform.cpp
  bsplines_non_uniform.cpp
  spline_1d.cpp
  spline_interpolator_1d.cpp
  spline_2d.cpp
  spline_interpolator_2d.cpp
)
target_include_directories(cxx_splines
    PUBLIC
        "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries( cxx_splines
  PUBLIC
     std::mdspan
     selalib::math_tools
     selalib::cxx_spline_matrix
     selalib::boundary_conditions
)

add_library(selalib::cxx_splines ALIAS cxx_splines)
