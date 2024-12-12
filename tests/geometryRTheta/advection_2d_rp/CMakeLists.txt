# SPDX-License-Identifier: MIT

#########################################################################################################
# All possibilities for MAPPING_TYPE:                                                                   #
# "CIRCULAR_MAPPING" "CZARNY_MAPPING_PHYSICAL" "CZARNY_MAPPING_PSEUDO_CARTESIAN" "DISCRETE_MAPPING"     #
# "REFINED_DISCRETE_MAPPING"                                                                            #
#                                                                                                       #
# All possiibilities for TIME_METHOD:                                                                   #
# "EULER_METHOD" "CRANK_NICOLSON_METHOD" "RK3_METHOD" "RK4_METHOD"                                      #
#                                                                                                       #
# All possibilities for SIMULATION:                                                                     #
# "TRANSLATION_SIMULATION" "ROTATION_SIMULATION" "DECENTRED_ROTATION_SIMULATION"                        #
#                                                                                                       #
# All possibilities for MESH: (NOT YET IMPLEMENTED)                                                     #
# "UNIFORM_MESH_R" "NON_UNIFORM_MESH_R_1" "NON_UNIFORM_MESH_R_2"                                        #
#                                                                                                       #
#                                                                                                       #
#########################################################################################################

set(TIME_METHOD "CRANK_NICOLSON_METHOD")

foreach(MAPPING_TYPE "CIRCULAR_MAPPING_PHYSICAL" "CZARNY_MAPPING_PHYSICAL" "CZARNY_MAPPING_PSEUDO_CARTESIAN" "DISCRETE_MAPPING_PSEUDO_CARTESIAN")
    foreach(SIMULATION "TRANSLATION_SIMULATION" "ROTATION_SIMULATION" "DECENTRED_ROTATION_SIMULATION")
        set(test_name "advection_${MAPPING_TYPE}__${TIME_METHOD}__${SIMULATION}")
        add_executable("${test_name}"
            test_cases.hpp
            advection_selected_test.cpp
        )
        target_link_libraries("${test_name}"
            PUBLIC
                DDC::core
                DDC::pdi
                paraconf::paraconf
                PDI::pdi
                Eigen3::Eigen

                sll::SLL

                gslx::advection_rp
                gslx::geometry_RTheta
                gslx::interpolation_2D_rp
                gslx::io
                gslx::paraconfpp
                gslx::quadrature
                gslx::timestepper
                gslx::utils
        )
        target_compile_definitions("${test_name}" PUBLIC -D${MAPPING_TYPE} -D${TIME_METHOD} -D${SIMULATION} -D${MESH})
        # CONVERGENCE ORDER TEST: ----------------------------------------
        add_test(NAME TestAdvectionRPConvergence_${MAPPING_TYPE}__${TIME_METHOD}__${SIMULATION}
               COMMAND "$<TARGET_FILE:Python3::Interpreter>" "${CMAKE_CURRENT_SOURCE_DIR}/display_all_errors_for_gtest.py"
                    "$<TARGET_FILE:advection_${MAPPING_TYPE}__${TIME_METHOD}__${SIMULATION}>")
        set_property(TEST TestAdvectionRPConvergence_${MAPPING_TYPE}__${TIME_METHOD}__${SIMULATION} PROPERTY TIMEOUT 50)
        set_property(TEST TestAdvectionRPConvergence_${MAPPING_TYPE}__${TIME_METHOD}__${SIMULATION} PROPERTY COST 10)
    endforeach()
endforeach()



##################################################################
# Create one executable for all the possibilties listed before.  #
#                                                                #
##################################################################

add_executable(advection_ALL
	test_cases.hpp
	advection_all_tests.cpp
)
target_link_libraries(advection_ALL
    PUBLIC
		DDC::core
		DDC::pdi
		paraconf::paraconf
		PDI::pdi
		Eigen3::Eigen

		sll::SLL

		gslx::advection_rp
		gslx::geometry_RTheta
		gslx::interpolation_2D_rp
		gslx::io
		gslx::paraconfpp
		gslx::quadrature
		gslx::timestepper
		gslx::utils
)


target_compile_definitions(advection_ALL PUBLIC)

find_package(Python3 REQUIRED COMPONENTS Interpreter)
