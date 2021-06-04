cmake_minimum_required(VERSION 3.15)

include(GoogleTest)

add_executable(ddc_tests
	block.cpp
        bsplines_non_uniform.cpp
        bsplines_uniform.cpp
        experimental/bsplines_uniform.cpp
        experimental/bsplines_non_uniform.cpp
        experimental/spline_builder.cpp
	mdomain.cpp
	mcoord.cpp
	nonuniformmesh.cpp
        nonuniformdomain.cpp
	taggedarray.cpp
	view.cpp
)
target_compile_features(ddc_tests PUBLIC cxx_std_17)
target_link_libraries(ddc_tests
	PUBLIC
		GTest::gtest_main
		ddc::ddc
)
gtest_discover_tests(ddc_tests)
