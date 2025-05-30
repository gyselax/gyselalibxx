


#########################################################################################################
# All possiibilities for SIMULATION:                                                                    #
# "TRANSLATION" "ROTATION" "DECENTERED_ROTATION"                                                        #
#                                                                                                       #
#########################################################################################################

function (advection_field_test_executable SIMULATION)
    set(test_name_gtest "advection_field_gtest_${SIMULATION}")

    add_executable("${test_name_gtest}"
            advection_field_gtest.cpp
            test_cases_adv_field.hpp
            ../../main.cpp
    )
    target_link_libraries("${test_name_gtest}"
    PUBLIC
            DDC::core
            DDC::pdi

            paraconf::paraconf

            PDI::pdi
            
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
    target_compile_definitions("${test_name_gtest}" PUBLIC -D${SIMULATION})    
    gtest_discover_tests("${test_name_gtest}"
        TEST_SUFFIX ${SIMULATION}
        DISCOVERY_MODE PRE_TEST)
endfunction()

# Select the parameters of the simulation: 
advection_field_test_executable(TRANSLATION)
advection_field_test_executable(ROTATION)
advection_field_test_executable(DECENTRED_ROTATION)

