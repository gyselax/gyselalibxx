cmake_minimum_required(VERSION 3.15)

include(GoogleTest)

add_executable(ddc_tests
	block.cpp
        bsplines_uniform.cpp
        bsplines_non_uniform.cpp
        spline_builder.cpp
	mdomain.cpp
	mcoord.cpp
        product_mdomain.cpp
        product_mesh.cpp
	non_uniform_mesh.cpp
	uniform_mesh.cpp
        single_mesh.cpp
	taggedarray.cpp
        type_seq.cpp
	view.cpp
)
target_compile_features(ddc_tests PUBLIC cxx_std_17)
target_link_libraries(ddc_tests
	PUBLIC
		GTest::gtest_main
		ddc::ddc
)
gtest_discover_tests(ddc_tests)
