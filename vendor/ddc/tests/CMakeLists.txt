cmake_minimum_required(VERSION 3.15)

include(GoogleTest)

add_executable(ddc_tests
	block.cpp
        deprecated/bsplines_non_uniform.cpp
        deprecated/bsplines_uniform.cpp
        bsplines_uniform.cpp
        bsplines_non_uniform.cpp
        spline_builder.cpp
        mesh_product.cpp
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
