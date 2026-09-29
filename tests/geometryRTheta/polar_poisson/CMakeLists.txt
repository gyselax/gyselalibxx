# SPDX-License-Identifier: MIT

foreach(MAPPING_TYPE "CIRCULAR_MAPPING" "CZARNY_MAPPING")
  foreach(SOLUTION "CURVILINEAR_SOLUTION" "CARTESIAN_SOLUTION")
    set(test_name "polar_poisson_convergence_${MAPPING_TYPE}_${SOLUTION}")
    add_executable("${test_name}"
        test_cases.cpp
        polarpoissonfemsolver.cpp
    )
    target_link_libraries("${test_name}"
        PUBLIC
            DDC::core
            DDC::pdi
            paraconf::paraconf
            PDI::PDI_C
            gyselalibxx::geometry_RTheta
            gyselalibxx::paraconfpp
            gyselalibxx::pde_solvers
            gyselalibxx::poisson_RTheta
            gyselalibxx::utils
    )
    target_compile_definitions("${test_name}" PUBLIC -D${MAPPING_TYPE} -D${SOLUTION})

    find_package(Python3 REQUIRED COMPONENTS Interpreter)

    if (NOT "${CMAKE_CXX_FLAGS}" MATCHES "--coverage")
        # Test is much too slow with coverage activated
        add_test(NAME TestPoissonConvergence_${MAPPING_TYPE}_${SOLUTION}
            COMMAND "$<TARGET_FILE:Python3::Interpreter>" "${CMAKE_CURRENT_SOURCE_DIR}/test_poisson.py"
                "$<TARGET_FILE:${test_name}>")
        set_property(TEST TestPoissonConvergence_${MAPPING_TYPE}_${SOLUTION} PROPERTY TIMEOUT 300)
        set_property(TEST TestPoissonConvergence_${MAPPING_TYPE}_${SOLUTION} PROPERTY COST 100)
    endif()
  endforeach()
endforeach()

add_executable(gmg_polar_poisson_like_solver_test
    gmgpolarpoissonsolver.cpp
    ../../main.cpp
)
target_link_libraries(gmg_polar_poisson_like_solver_test
    PUBLIC
        GTest::gtest
	    GTest::gmock
        DDC::core
        DDC::splines
        gyselalibxx::geometry_RTheta
        gyselalibxx::pde_solvers
)

gtest_discover_tests(gmg_polar_poisson_like_solver_test
    DISCOVERY_TIMEOUT ${GYSELALIBXX_DISCOVERY_TIMEOUT}
    DISCOVERY_MODE PRE_TEST)
