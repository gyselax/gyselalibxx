cmake_minimum_required(VERSION 3.15)

include(GoogleTest)

set(SLL_SPLINES_TEST_DEGREE_MIN 3 CACHE STRING "Minimum degree to test splines.")
set(SLL_SPLINES_TEST_DEGREE_MAX 3 CACHE STRING "Maximum degree to test splines.")

add_executable(misc_tests_sll
    main.cpp
    view.cpp
    gauss_legendre_integration.cpp
)
target_compile_features(misc_tests_sll PUBLIC cxx_std_17)
target_link_libraries(misc_tests_sll
    PUBLIC
        GTest::gtest
        sll::SLL
)

add_executable(bsplines_tests_sll
	main.cpp
	polar_bsplines.cpp
)
target_compile_features(bsplines_tests_sll PUBLIC cxx_std_17)
target_link_libraries(bsplines_tests_sll
    PUBLIC
        GTest::gtest
        sll::SLL
)

gtest_discover_tests(bsplines_tests_sll)

add_executable(matrix_tests_sll
  main.cpp
  matrix.cpp
  matrix_batch_ell.cpp
  matrix_batch_tridiag.cpp
)
target_compile_features(matrix_tests_sll PUBLIC cxx_std_17)
target_link_libraries(matrix_tests_sll
    PUBLIC
        GTest::gtest
        sll::SLL
)
gtest_discover_tests(matrix_tests_sll)

add_executable(bernstein_tests_sll
  main.cpp
	bernstein.cpp
)
target_compile_features(bernstein_tests_sll PUBLIC cxx_std_17)
target_link_libraries(bernstein_tests_sll
    PUBLIC
        GTest::gtest
        sll::SLL
)
gtest_discover_tests(bernstein_tests_sll)

add_executable(metric_tensor_tests_sll
  main.cpp
	metric_tensor.cpp
)
target_compile_features(metric_tensor_tests_sll PUBLIC cxx_std_17)
target_link_libraries(metric_tensor_tests_sll
    PUBLIC
        GTest::gtest
        sll::SLL
)
gtest_discover_tests(metric_tensor_tests_sll)

add_executable(jacobian_mapping_tests_sll
  main.cpp
	mapping_jacobian.cpp
)
target_compile_features(jacobian_mapping_tests_sll PUBLIC cxx_std_17)
target_link_libraries(jacobian_mapping_tests_sll
    PUBLIC
        GTest::gtest
        sll::SLL
)
gtest_discover_tests(jacobian_mapping_tests_sll)

add_executable(jacobian_mapping_matrix_coef_tests_sll
  main.cpp
	mapping_jacobian_matrix_coef.cpp
)
target_compile_features(jacobian_mapping_matrix_coef_tests_sll PUBLIC cxx_std_17)
target_link_libraries(jacobian_mapping_matrix_coef_tests_sll
    PUBLIC
        GTest::gtest
        sll::SLL
)
gtest_discover_tests(jacobian_mapping_matrix_coef_tests_sll)

add_executable(pseudo_cartesian_tests_sll
  main.cpp
	pseudo_cartesian_matrix.cpp
)
target_compile_features(pseudo_cartesian_tests_sll PUBLIC cxx_std_17)
target_link_libraries(pseudo_cartesian_tests_sll
    PUBLIC
        GTest::gtest
        sll::SLL
)
gtest_discover_tests(pseudo_cartesian_tests_sll)

add_executable(refined_discrete_mapping_test
  main.cpp
	refined_discrete_mapping.cpp
)
target_compile_features(refined_discrete_mapping_test PUBLIC cxx_std_17)
target_link_libraries(refined_discrete_mapping_test
    PUBLIC
        GTest::gtest
        sll::SLL
)
gtest_discover_tests(refined_discrete_mapping_test)

foreach(CONTINUITY RANGE -1 1)
  math(EXPR MIN_DEGREE "${CONTINUITY}+1")
  if (${MIN_DEGREE} LESS 1)
    set(MIN_DEGREE 1)
  endif()
  if (${MIN_DEGREE} LESS ${SLL_SPLINES_TEST_DEGREE_MIN})
    set(MIN_DEGREE ${SLL_SPLINES_TEST_DEGREE_MIN})
  endif()
  foreach(DEGREE_R RANGE ${MIN_DEGREE} ${SLL_SPLINES_TEST_DEGREE_MAX})
    foreach(DEGREE_P RANGE ${MIN_DEGREE} ${SLL_SPLINES_TEST_DEGREE_MAX})
      foreach(BSPLINES_TYPE "BSPLINES_TYPE_UNIFORM" "BSPLINES_TYPE_NON_UNIFORM")
        foreach(MAPPING_TYPE "CIRCULAR_MAPPING" "CZARNY_MAPPING")
          set(test_name "polar_splines_tests_sll_DEGREE_R_${DEGREE_R}_DEGREE_P_${DEGREE_P}_CONTINUITY_${CONTINUITY}_${BSPLINES_TYPE}_${MAPPING_TYPE}")
          add_executable("${test_name}" polar_splines.cpp)
          target_compile_features("${test_name}" PUBLIC cxx_std_17)
          target_link_libraries("${test_name}"
            PUBLIC
            GTest::gtest
            sll::SLL
          )
          target_compile_definitions("${test_name}" PUBLIC -DDEGREE_R=${DEGREE_R} -DDEGREE_P=${DEGREE_P} -D${BSPLINES_TYPE} -DCONTINUITY=${CONTINUITY} -D${MAPPING_TYPE})
          add_test("${test_name}" "${test_name}")
        endforeach()
      endforeach()
    endforeach()
  endforeach()
endforeach()
