

#########################################################################################################
# All possiibilities for TIME_METHOD:                                                                   #
# "EULER_METHOD" "CRANK_NICOLSON_METHOD" "RK3_METHOD" "RK4_METHOD"                                      #
#                                                                                                       #
# All possibilities for PREDCORR:                                                                       #
# "PREDCORR" "EXPLICIT_PREDCORR" "IMPLICIT_PREDCORR"                                                    #
#                                                                                                       #
#########################################################################################################


function (diocotron_executable PREDCORR_METHOD TIME_METHOD )
    set(test_name "diocotron_${PREDCORR_METHOD}_${TIME_METHOD}")

    add_executable("${test_name}"
         diocotron.cpp
    )
    target_link_libraries("${test_name}"
        PUBLIC
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
    target_compile_definitions("${test_name}" PUBLIC -D${TIME_METHOD} -D${PREDCORR_METHOD})
endfunction()


# Select the parameters of the simulation:
diocotron_executable(PREDCORR EULER_METHOD)
diocotron_executable(PREDCORR CRANK_NICOLSON_METHOD)
diocotron_executable(PREDCORR RK3_METHOD)
diocotron_executable(PREDCORR RK4_METHOD)

diocotron_executable(EXPLICIT_PREDCORR EULER_METHOD)
diocotron_executable(IMPLICIT_PREDCORR EULER_METHOD)


