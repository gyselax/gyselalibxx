cmake_minimum_required(VERSION 3.15)


add_executable(Li_norms_spline_quadrature_tests
	../../main.cpp
	tests_L1_and_L2_norms.cpp
)
target_compile_features(Li_norms_spline_quadrature_tests PUBLIC cxx_std_17)
target_link_libraries(Li_norms_spline_quadrature_tests
    PUBLIC
        GTest::gtest
        GTest::gmock
        DDC::DDC
        sll::SLL
        gslx::geometry_RTheta
        gslx::quadrature
)
gtest_discover_tests(Li_norms_spline_quadrature_tests)
