


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
            DDC::DDC
            DDC::PDI_Wrapper
            paraconf::paraconf
            PDI::pdi
            
            Eigen3::Eigen
            
            sll::SLL
            gslx::paraconfpp

            GTest::gtest
            GTest::gmock
            
            gslx::poisson_RTheta
            gslx::geometry_RTheta
            gslx::interpolation_2D_rp
            gslx::advection_rp
            gslx::time_integration_rp
            gslx::quadrature
            gslx::utils
            gslx::diocotron_init_eq
            gslx::simulation_utils
    )
    target_compile_definitions("${test_name_gtest}" PUBLIC -D${SIMULATION})    
    gtest_discover_tests("${test_name_gtest}" DISCOVERY_MODE PRE_TEST)
endfunction()

# Select the parameters of the simulation: 
advection_field_test_executable(TRANSLATION)
advection_field_test_executable(ROTATION)
advection_field_test_executable(DECENTRED_ROTATION)

