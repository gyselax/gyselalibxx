
foreach(PREDCORR "RK2_PREDCORR" "EXPLICIT_PREDCORR" "IMPLICIT_PREDCORR")
    set(test_name "bsl_${PREDCORR}_test")
    add_executable("${test_name}"
        ../../main.cpp
        bsl_predcorr_test.cpp
    )
    target_link_libraries("${test_name}"
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
    target_compile_definitions("${test_name}" PUBLIC -D${PREDCORR} )
    gtest_discover_tests("${test_name}" DISCOVERY_MODE PRE_TEST)
endforeach()