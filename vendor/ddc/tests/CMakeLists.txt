# SPDX-License-Identifier: MIT
cmake_minimum_required(VERSION 3.15)

include(GoogleTest)

add_executable(ddc_tests
	block.cpp
	mdomain.cpp
	mcoord.cpp
	product_mdomain.cpp
	product_mesh.cpp
	non_uniform_mesh.cpp
	uniform_mesh.cpp
	single_mesh.cpp
	taggedvector.cpp
	type_seq.cpp
)
target_compile_features(ddc_tests PUBLIC cxx_std_17)
target_link_libraries(ddc_tests
	PUBLIC
		GTest::gtest_main
		DDC::DDC
)
gtest_discover_tests(ddc_tests)
