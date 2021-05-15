cmake_minimum_required(VERSION 3.15)

include(GoogleTest)

add_executable(unit_tests
	blockview.cpp
	mdomain.cpp
	mcoord.cpp
	taggedarray.cpp
	splines.cpp
	view.cpp
    fft.cpp
)
target_compile_features(unit_tests PUBLIC cxx_std_17)
target_link_libraries(unit_tests
	PUBLIC
		GTest::gtest_main
		sll::splines
		vcx::array_tools
		vcx::mesh
		FFTW::Double
)
gtest_discover_tests(unit_tests)
