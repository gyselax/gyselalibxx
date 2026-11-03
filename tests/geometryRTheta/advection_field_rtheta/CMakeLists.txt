


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
        
        gyselalibxx::advection
        gyselalibxx::diocotron_init_eq
        gyselalibxx::geometry_RTheta
        gyselalibxx::interpolation
        gyselalibxx::paraconfpp
        gyselalibxx::poisson_RTheta
        gyselalibxx::quadrature
        gyselalibxx::simulation_utils
        gyselalibxx::time_integration_RTheta
        gyselalibxx::utils
)
gtest_discover_tests("${test_name_gtest}"
    DISCOVERY_TIMEOUT ${GYSELALIBXX_DISCOVERY_TIMEOUT}
    DISCOVERY_MODE PRE_TEST)

