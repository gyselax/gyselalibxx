cmake_minimum_required(VERSION 3.15)

include(GoogleTest)

set(SLL_SPLINES_TEST_DEGREE_MIN 3 CACHE STRING "Minimum degree to test splines.")
set(SLL_SPLINES_TEST_DEGREE_MAX 3 CACHE STRING "Maximum degree to test splines.")

add_executable(splines_tests_sll
    main.cpp
    chunk.cpp
    view.cpp
    gauss_legendre_integration.cpp
)
if(GYSELALIBXX_ENABLE_DEPRECATED)
    target_sources(splines_tests_sll
	PRIVATE
        bsplines_non_uniform.cpp
        bsplines_uniform.cpp
        splines.cpp
    )
endif()
target_compile_features(splines_tests_sll PUBLIC cxx_std_17)
target_link_libraries(splines_tests_sll
    PUBLIC
        GTest::gtest
        sll::splines
)
gtest_discover_tests(splines_tests_sll)

add_executable(bsplines_tests_sll
	main.cpp
	bsplines.cpp
	polar_bsplines.cpp
)
target_compile_features(bsplines_tests_sll PUBLIC cxx_std_17)
target_link_libraries(bsplines_tests_sll
    PUBLIC
        GTest::gtest
        sll::splines
)

gtest_discover_tests(bsplines_tests_sll)

add_executable(matrix_tests_sll
  main.cpp
	matrix.cpp
)
target_compile_features(matrix_tests_sll PUBLIC cxx_std_17)
target_link_libraries(matrix_tests_sll
    PUBLIC
        GTest::gtest
        sll::splines
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
        sll::splines
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
        sll::splines
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
        sll::splines
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
        sll::splines
)
gtest_discover_tests(jacobian_mapping_matrix_coef_tests_sll)


add_executable(const_extrapol_tests_sll
  	main.cpp
	constant_extrapolation_bc.cpp
)
target_compile_features(const_extrapol_tests_sll PUBLIC cxx_std_17)
target_link_libraries(const_extrapol_tests_sll
    PUBLIC
        GTest::gtest
        sll::splines
)
gtest_discover_tests(const_extrapol_tests_sll)



add_executable(pseudo_cartesian_tests_sll
  	main.cpp
	pseudo_cartesian_matrix.cpp
)
target_compile_features(pseudo_cartesian_tests_sll PUBLIC cxx_std_17)
target_link_libraries(pseudo_cartesian_tests_sll
    PUBLIC
        GTest::gtest
        sll::splines
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
        sll::splines
)
gtest_discover_tests(refined_discrete_mapping_test)




foreach(DEGREE_X RANGE "${SLL_SPLINES_TEST_DEGREE_MIN}" "${SLL_SPLINES_TEST_DEGREE_MAX}")
  foreach(BSPLINES_TYPE "BSPLINES_TYPE_UNIFORM" "BSPLINES_TYPE_NON_UNIFORM")
    set(test_name "splines_tests_sll_DEGREE_X_${DEGREE_X}_${BSPLINES_TYPE}")
    add_executable("${test_name}" periodic_spline_builder.cpp)
    target_compile_features("${test_name}" PUBLIC cxx_std_17)
    target_link_libraries("${test_name}"
      PUBLIC
      GTest::gtest
      sll::splines
    )
    target_compile_definitions("${test_name}" PUBLIC -DDEGREE_X=${DEGREE_X} -D${BSPLINES_TYPE})
    add_test("${test_name}" "${test_name}")
  endforeach()
endforeach()

foreach(DEGREE_X RANGE "${SLL_SPLINES_TEST_DEGREE_MIN}" "${SLL_SPLINES_TEST_DEGREE_MAX}")
  set(test_name "splines_ordered_points_tests_sll_DEGREE_X_${DEGREE_X}")
  add_executable("${test_name}" periodic_spline_builder_ordered_points.cpp)
  target_compile_features("${test_name}" PUBLIC cxx_std_17)
  target_link_libraries("${test_name}"
    PUBLIC
    GTest::gtest
    sll::splines
  )
  target_compile_definitions("${test_name}" PUBLIC -DDEGREE_X=${DEGREE_X})
  add_test("${test_name}" "${test_name}")
endforeach()

foreach(BCL "BCL_GREVILLE" "BCL_HERMITE")
  foreach(BCR "BCR_GREVILLE" "BCR_HERMITE")
    foreach(EVALUATOR "EVALUATOR_COSINE" "EVALUATOR_POLYNOMIAL")
      foreach(DEGREE_X RANGE "${SLL_SPLINES_TEST_DEGREE_MIN}" "${SLL_SPLINES_TEST_DEGREE_MAX}")
        foreach(BSPLINES_TYPE "BSPLINES_TYPE_UNIFORM" "BSPLINES_TYPE_NON_UNIFORM")
          set(test_name "splines_tests_sll_DEGREE_X_${DEGREE_X}_${BSPLINES_TYPE}_${EVALUATOR}_bcl_${BCL}_bcr_${BCR}")
          add_executable("${test_name}" non_periodic_spline_builder.cpp)
          target_compile_features("${test_name}" PUBLIC cxx_std_17)
          target_link_libraries("${test_name}"
            PUBLIC
            GTest::gtest
            sll::splines
          )
          target_compile_definitions("${test_name}" PUBLIC -DDEGREE_X=${DEGREE_X} -D${BSPLINES_TYPE} -D${EVALUATOR} -D${BCL} -D${BCR})
          add_test("${test_name}" "${test_name}")
        endforeach()
      endforeach()
    endforeach()
  endforeach()
endforeach()

foreach(DEGREE_X RANGE "${SLL_SPLINES_TEST_DEGREE_MIN}" "${SLL_SPLINES_TEST_DEGREE_MAX}")
  foreach(DEGREE_Y RANGE "${SLL_SPLINES_TEST_DEGREE_MIN}" "${SLL_SPLINES_TEST_DEGREE_MAX}")
    foreach(BSPLINES_TYPE "BSPLINES_TYPE_UNIFORM" "BSPLINES_TYPE_NON_UNIFORM")
      set(test_name "2d_splines_tests_sll_DEGREE_X_${DEGREE_X}_DEGREE_Y_${DEGREE_Y}_${BSPLINES_TYPE}")
      add_executable("${test_name}" 2d_spline_builder_periodic.cpp)
      target_compile_features("${test_name}" PUBLIC cxx_std_17)
      target_link_libraries("${test_name}"
        PUBLIC
        GTest::gtest
        sll::splines
      )
      target_compile_definitions("${test_name}" PUBLIC -DDEGREE_X=${DEGREE_X} -DDEGREE_Y=${DEGREE_Y} -D${BSPLINES_TYPE})
      add_test("${test_name}" "${test_name}")
    endforeach()
  endforeach()
endforeach()

foreach(BCL "BCL_GREVILLE" "BCL_HERMITE")
  foreach(BCR "BCR_GREVILLE" "BCR_HERMITE")
    foreach(EVALUATOR "EVALUATOR_POLYNOMIAL")
      foreach(DEGREE_X RANGE "${SLL_SPLINES_TEST_DEGREE_MIN}" "${SLL_SPLINES_TEST_DEGREE_MAX}")
        foreach(DEGREE_Y RANGE "${SLL_SPLINES_TEST_DEGREE_MIN}" "${SLL_SPLINES_TEST_DEGREE_MAX}")
          foreach(BSPLINES_TYPE "BSPLINES_TYPE_UNIFORM" "BSPLINES_TYPE_NON_UNIFORM")
            set(test_name "2d_splines_tests_sll_DEGREE_X_${DEGREE_X}_DEGREE_Y_${DEGREE_Y}_${BSPLINES_TYPE}_${EVALUATOR}_${BCL}_${BCR}")
            add_executable("${test_name}" 2d_spline_builder.cpp)
            target_compile_features("${test_name}" PUBLIC cxx_std_17)
            target_link_libraries("${test_name}"
              PUBLIC
              GTest::gtest
              sll::splines
            )
            target_compile_definitions("${test_name}" PUBLIC -DDEGREE_X=${DEGREE_X} -DDEGREE_Y=${DEGREE_Y} -D${BSPLINES_TYPE} -D${EVALUATOR} -D${BCL} -D${BCR})
            add_test("${test_name}" "${test_name}")
          endforeach()
        endforeach()
      endforeach()
    endforeach()
  endforeach()
endforeach()

foreach(CONTINUITY RANGE -1 1)
  math(EXPR MIN_DEGREE "${CONTINUITY}+1")
  if (${MIN_DEGREE} LESS 1)
    set(MIN_DEGREE 1)
  endif()
  foreach(DEGREE_R RANGE ${MIN_DEGREE} 6)
    foreach(DEGREE_P RANGE ${MIN_DEGREE} 6)
      foreach(BSPLINES_TYPE "BSPLINES_TYPE_UNIFORM" "BSPLINES_TYPE_NON_UNIFORM")
        foreach(MAPPING_TYPE "CIRCULAR_MAPPING" "CZARNY_MAPPING")
          set(test_name "polar_splines_tests_sll_DEGREE_R_${DEGREE_R}_DEGREE_P_${DEGREE_P}_CONTINUITY_${CONTINUITY}_${BSPLINES_TYPE}_${MAPPING_TYPE}")
          add_executable("${test_name}" polar_splines.cpp)
          target_compile_features("${test_name}" PUBLIC cxx_std_17)
          target_link_libraries("${test_name}"
            PUBLIC
            GTest::gtest
            sll::splines
          )
          target_compile_definitions("${test_name}" PUBLIC -DDEGREE_R=${DEGREE_R} -DDEGREE_P=${DEGREE_P} -D${BSPLINES_TYPE} -DCONTINUITY=${CONTINUITY} -D${MAPPING_TYPE})
          add_test("${test_name}" "${test_name}")
        endforeach()
      endforeach()
    endforeach()
  endforeach()
endforeach()
