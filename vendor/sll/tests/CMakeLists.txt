cmake_minimum_required(VERSION 3.15)

include(GoogleTest)

add_executable(splines_tests
    main.cpp
    chunk.cpp
    view.cpp
    gauss_legendre_integration.cpp
)
if(VOICEXX_ENABLE_DEPRECATED)
    target_sources(splines_tests
	PRIVATE
        bsplines_non_uniform.cpp
        bsplines_uniform.cpp
        splines.cpp
    )
endif()
target_compile_features(splines_tests PUBLIC cxx_std_17)
target_link_libraries(splines_tests
    PUBLIC
        GTest::gtest
        sll::splines
)
gtest_discover_tests(splines_tests)

add_executable(bsplines_tests
  main.cpp
	bsplines.cpp
)
target_compile_features(bsplines_tests PUBLIC cxx_std_17)
target_link_libraries(bsplines_tests
    PUBLIC
        GTest::gtest
        sll::splines
)

gtest_discover_tests(bsplines_tests)

add_executable(matrix_tests
  main.cpp
	matrix.cpp
)
target_compile_features(matrix_tests PUBLIC cxx_std_17)
target_link_libraries(matrix_tests
    PUBLIC
        GTest::gtest
        sll::splines
)
gtest_discover_tests(matrix_tests)

add_executable(bernstein_tests
  main.cpp
	bernstein.cpp
)
target_compile_features(bernstein_tests PUBLIC cxx_std_17)
target_link_libraries(bernstein_tests
    PUBLIC
        GTest::gtest
        sll::splines
)
gtest_discover_tests(bernstein_tests)

foreach(EVALUATOR "Cosine" "Polynomial")
  foreach(DEGREE_X RANGE 1 7)
    foreach(UNIFORM 0 1)
      set(test_name "splines_tests_degree_${DEGREE_X}_uniform_${UNIFORM}_evaluator_${EVALUATOR}")
      add_executable("${test_name}" periodic_spline_builder.cpp)
      target_compile_features("${test_name}" PUBLIC cxx_std_17)
      target_link_libraries("${test_name}"
        PUBLIC
        GTest::gtest
        sll::splines
      )
      target_compile_definitions("${test_name}" PUBLIC -DDEGREE_X=${DEGREE_X} -DUNIFORM=${UNIFORM} -DEVALUATOR=${EVALUATOR})
      add_test("${test_name}" "${test_name}")
    endforeach()
  endforeach()
endforeach()

foreach(DEGREE_X RANGE 1 7)
  set(test_name "splines_ordered_points_degree_${DEGREE_X}")
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

foreach(BCL "GREVILLE" "HERMITE")
  foreach(BCR "GREVILLE" "HERMITE")
    foreach(EVALUATOR "Cosine" "Polynomial")
      foreach(DEGREE_X RANGE 1 7)
        foreach(UNIFORM 0 1)
          set(test_name "splines_tests_degree_${DEGREE_X}_uniform_${UNIFORM}_evaluator_${EVALUATOR}_bcl_${BCL}_bcr_${BCR}")
          add_executable("${test_name}" non_periodic_spline_builder.cpp)
          target_compile_features("${test_name}" PUBLIC cxx_std_17)
          target_link_libraries("${test_name}"
            PUBLIC
            GTest::gtest
            sll::splines
          )
          target_compile_definitions("${test_name}" PUBLIC -DDEGREE_X=${DEGREE_X} -DUNIFORM=${UNIFORM} -DEVALUATOR=${EVALUATOR} -DBCL=${BCL} -DBCR=${BCR})
          add_test("${test_name}" "${test_name}")
        endforeach()
      endforeach()
    endforeach()
  endforeach()
endforeach()

foreach(EVALUATOR "Cosine")
  foreach(DEGREE_X RANGE 1 7)
    foreach(DEGREE_Y RANGE 1 7)
      foreach(UNIFORM 0 1)
        set(test_name "2d_splines_tests_degree_x_${DEGREE_X}_degree_y_${DEGREE_Y}_uniform_${UNIFORM}_evaluator_${EVALUATOR}")
        add_executable("${test_name}" 2d_spline_builder_periodic.cpp)
        target_compile_features("${test_name}" PUBLIC cxx_std_17)
        target_link_libraries("${test_name}"
          PUBLIC
          GTest::gtest
          sll::splines
        )
        target_compile_definitions("${test_name}" PUBLIC -DDEGREE_X=${DEGREE_X} -DDEGREE_Y=${DEGREE_Y} -DUNIFORM=${UNIFORM} -DEVALUATOR=${EVALUATOR})
        add_test("${test_name}" "${test_name}")
      endforeach()
    endforeach()
  endforeach()
endforeach()

foreach(BCL "GREVILLE" "HERMITE")
  foreach(BCR "GREVILLE" "HERMITE")
    foreach(EVALUATOR "Polynomial")
      foreach(DEGREE_X RANGE 1 7)
        foreach(DEGREE_Y RANGE 1 7)
          foreach(UNIFORM 0 1)
            set(test_name "2d_splines_tests_degree_x_${DEGREE_X}_degree_y_${DEGREE_Y}_uniform_${UNIFORM}_evaluator_${EVALUATOR}_bcl_${BCL}_bcr_${BCR}")
            add_executable("${test_name}" 2d_spline_builder.cpp)
            target_compile_features("${test_name}" PUBLIC cxx_std_17)
            target_link_libraries("${test_name}"
              PUBLIC
              GTest::gtest
              sll::splines
            )
            target_compile_definitions("${test_name}" PUBLIC -DDEGREE_X=${DEGREE_X} -DDEGREE_Y=${DEGREE_Y} -DUNIFORM=${UNIFORM} -DEVALUATOR=${EVALUATOR} -DBCL=${BCL} -DBCR=${BCR})
            add_test("${test_name}" "${test_name}")
          endforeach()
        endforeach()
      endforeach()
    endforeach()
  endforeach()
endforeach()
