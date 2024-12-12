

add_executable(Li_norms_spline_quadrature_tests
	../../main.cpp
	tests_L1_and_L2_norms.cpp
)
target_link_libraries(Li_norms_spline_quadrature_tests
    PUBLIC
        GTest::gtest
        GTest::gmock
        DDC::DDC
        sll::SLL
        gslx::geometry_RTheta
        gslx::quadrature
        gslx::math_tools
        gslx::utils

)
gtest_discover_tests(Li_norms_spline_quadrature_tests DISCOVERY_MODE PRE_TEST)
