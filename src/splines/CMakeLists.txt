add_library(splines STATIC
	bsplines.cpp
	bsplines_uniform.cpp
	bsplines_non_uniform.cpp
	spline_1d.cpp
	spline_interpolator_1d.cpp
	spline_2d.cpp
	spline_interpolator_2d.cpp
)

target_compile_features(splines
	PUBLIC
		cxx_std_17
)

target_include_directories(splines
	PUBLIC
		"${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries(splines
	PUBLIC
		std::mdspan
		sll::math_tools
		sll::matrix
		sll::boundary_conditions
)

add_library(sll::splines ALIAS splines)
