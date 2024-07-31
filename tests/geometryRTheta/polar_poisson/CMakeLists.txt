
foreach(MAPPING_TYPE "CIRCULAR_MAPPING" "CZARNY_MAPPING")
  foreach(SOLUTION "CURVILINEAR_SOLUTION" "CARTESIAN_SOLUTION")
    set(test_name "polar_poisson_convergence_${MAPPING_TYPE}_${SOLUTION}")
    add_executable("${test_name}"
        test_cases.cpp
        polarpoissonfemsolver.cpp
    )
    target_link_libraries("${test_name}"
        PUBLIC
        DDC::DDC
            DDC::PDI_Wrapper
            paraconf::paraconf
            PDI::pdi
            sll::SLL
            gslx::poisson_RTheta
            gslx::paraconfpp
            Eigen3::Eigen
            gslx::geometry_RTheta
            gslx::utils

    )
    target_compile_definitions("${test_name}" PUBLIC -D${MAPPING_TYPE} -D${SOLUTION})

    find_package(Python3 REQUIRED COMPONENTS Interpreter)

    add_test(NAME TestPoissonConvergence_${MAPPING_TYPE}_${SOLUTION}
        COMMAND "$<TARGET_FILE:Python3::Interpreter>" "${CMAKE_CURRENT_SOURCE_DIR}/test_poisson.py"
            "$<TARGET_FILE:${test_name}>")
    set_property(TEST TestPoissonConvergence_${MAPPING_TYPE}_${SOLUTION} PROPERTY TIMEOUT 200)
    set_property(TEST TestPoissonConvergence_${MAPPING_TYPE}_${SOLUTION} PROPERTY COST 100)
  endforeach()
endforeach()
