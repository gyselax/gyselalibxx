

add_executable(predcorr_test
	../../main.cpp
	bsl_predcorr_second_order_explicit_test.cpp
	# bsl_predcorr_second_order_implicit_test.cpp
	# bsl_predcorr_test.cpp
)
target_link_libraries(predcorr_test
    PUBLIC
        GTest::gtest
        GTest::gmock
        DDC::core
        DDC::pdi

        paraconf::paraconf

        PDI::pdi

        gslx::advection
        gslx::diocotron_init_eq
        gslx::geometry_RTheta
        gslx::interpolation
        gslx::io
        gslx::paraconfpp
        gslx::poisson_RTheta
        gslx::quadrature
        gslx::simulation_utils
        gslx::time_integration_RTheta
        gslx::utils
)
gtest_discover_tests(predcorr_test DISCOVERY_MODE PRE_TEST)
