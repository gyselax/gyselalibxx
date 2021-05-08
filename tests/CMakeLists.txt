cmake_minimum_required(VERSION 3.15)

include(GoogleTest)

add_executable(unit_tests
	block_test.cpp
	mdomain_test.cpp
	taggedarray_test.cpp
)
target_compile_features(unit_tests PUBLIC cxx_std_17)
target_link_libraries(unit_tests
	PUBLIC
		GTest::gtest_main
		vcx::array_tools
		vcx::mesh
)
gtest_discover_tests(unit_tests)
