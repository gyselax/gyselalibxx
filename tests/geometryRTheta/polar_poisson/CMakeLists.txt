
foreach(MAPPING_TYPE "CIRCULAR_MAPPING" "CZARNY_MAPPING")
  foreach(SOLUTION "CURVILINEAR_SOLUTION" "CARTESIAN_SOLUTION")
    set(test_name "polar_poisson_convergence_${MAPPING_TYPE}_${SOLUTION}")
    add_executable("${test_name}"
        test_cases.cpp
        polarpoissonfemsolver.cpp
    )
    target_compile_features("${test_name}" PUBLIC cxx_std_17)
    target_link_libraries("${test_name}"
        PUBLIC
        DDC::DDC
            DDC::PDI_Wrapper
            paraconf::paraconf
            PDI::pdi
            sll::splines
            gslx::poisson_RTheta
            gslx::paraconfpp
            Eigen3::Eigen
            gslx::geometry_RTheta
    )
    target_compile_definitions("${test_name}" PUBLIC -D${MAPPING_TYPE} -D${SOLUTION})

    find_package(Python3 REQUIRED COMPONENTS Interpreter)

    add_test(NAME TestPoissonConvergence_${MAPPING_TYPE}_${SOLUTION}
        COMMAND "$<TARGET_FILE:Python3::Interpreter>" "${CMAKE_CURRENT_SOURCE_DIR}/test_poisson.py"
            "$<TARGET_FILE:${test_name}>")
    set_property(TEST TestPoissonConvergence_${MAPPING_TYPE}_${SOLUTION} PROPERTY TIMEOUT 200)
  endforeach()
endforeach()




foreach(MAPPING_TYPE "CIRCULAR_MAPPING" "CZARNY_MAPPING")
    set(test_name "vlasov_polar_poisson_convergence_${MAPPING_TYPE}")
    add_executable("${test_name}"
        test_cases.cpp
        vlasovpoissonsolver.cpp
    )
    target_compile_features("${test_name}" PUBLIC cxx_std_17)
    target_link_libraries("${test_name}"
        PUBLIC
        DDC::DDC
            DDC::PDI_Wrapper
            paraconf::paraconf
            PDI::pdi
            sll::splines
            gslx::poisson_RTheta
            gslx::paraconfpp
            Eigen3::Eigen
            gslx::geometry_RTheta
    )
    target_compile_definitions("${test_name}" PUBLIC -D${MAPPING_TYPE})

    find_package(Python3 REQUIRED COMPONENTS Interpreter)

    add_test(NAME TestPoissonConvergence_${MAPPING_TYPE}_${SOLUTION}
        COMMAND "$<TARGET_FILE:Python3::Interpreter>" "${CMAKE_CURRENT_SOURCE_DIR}/test_vlasov_poisson.py"
            "$<TARGET_FILE:${test_name}>")
    set_property(TEST TestPoissonConvergence_${MAPPING_TYPE}_${SOLUTION} PROPERTY TIMEOUT 200)
endforeach()
