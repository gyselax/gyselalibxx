add_library(splines STATIC
	spline_builder_1d.cpp
	spline_builder_2d.cpp
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
		ddc::ddc
		sll::matrix
)

add_library(sll::splines ALIAS splines)
