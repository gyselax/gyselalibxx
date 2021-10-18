# SPDX-License-Identifier: MIT
cmake_minimum_required(VERSION 3.15)

include(GoogleTest)

add_executable(ddc_tests
	chunck.cpp
	discrete_domain.cpp
	discrete_space.cpp
	non_uniform_discretization.cpp
	single_discretization.cpp
	tagged_vector.cpp
	type_seq.cpp
	uniform_discretization.cpp
)
target_compile_features(ddc_tests PUBLIC cxx_std_17)
target_link_libraries(ddc_tests
	PUBLIC
		GTest::gtest_main
		DDC::DDC
)
gtest_discover_tests(ddc_tests)
