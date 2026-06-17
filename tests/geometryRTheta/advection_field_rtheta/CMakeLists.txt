


set(test_name_gtest "advection_field_gtest")

add_executable("${test_name_gtest}"
        advection_field_gtest.cpp
        ../../main.cpp
)
target_link_libraries("${test_name_gtest}"
PUBLIC
        DDC::core
        DDC::pdi

        paraconf::paraconf

        PDI::PDI_C
        
        GTest::gtest
        GTest::gmock
        
        gslx::advection
        gslx::diocotron_init_eq
        gslx::geometry_RTheta
        gslx::interpolation
        gslx::paraconfpp
        gslx::poisson_RTheta
        gslx::quadrature
        gslx::simulation_utils
        gslx::time_integration_RTheta
        gslx::utils
)
gtest_discover_tests("${test_name_gtest}"
    DISCOVERY_TIMEOUT ${GYSELALIBXX_DISCOVERY_TIMEOUT}
    DISCOVERY_MODE PRE_TEST)

