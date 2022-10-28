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

foreach(DEGREE_X RANGE 1 7)
  foreach(BSPLINES_TYPE "BSPLINES_TYPE_UNIFORM" "BSPLINES_TYPE_NON_UNIFORM")
    set(test_name "splines_tests_DEGREE_X_${DEGREE_X}_${BSPLINES_TYPE}")
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

foreach(DEGREE_X RANGE 1 7)
  set(test_name "splines_ordered_points_DEGREE_X_${DEGREE_X}")
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
      foreach(DEGREE_X RANGE 1 7)
        foreach(BSPLINES_TYPE "BSPLINES_TYPE_UNIFORM" "BSPLINES_TYPE_NON_UNIFORM")
          set(test_name "splines_tests_DEGREE_X_${DEGREE_X}_${BSPLINES_TYPE}_${EVALUATOR}_bcl_${BCL}_bcr_${BCR}")
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

foreach(DEGREE_X RANGE 1 7)
  foreach(DEGREE_Y RANGE 1 7)
    foreach(BSPLINES_TYPE "BSPLINES_TYPE_UNIFORM" "BSPLINES_TYPE_NON_UNIFORM")
      set(test_name "2d_splines_tests_DEGREE_X_${DEGREE_X}_DEGREE_Y_${DEGREE_Y}_${BSPLINES_TYPE}")
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
      foreach(DEGREE_X RANGE 1 7)
        foreach(DEGREE_Y RANGE 1 7)
          foreach(BSPLINES_TYPE "BSPLINES_TYPE_UNIFORM" "BSPLINES_TYPE_NON_UNIFORM")
            set(test_name "2d_splines_tests_DEGREE_X_${DEGREE_X}_DEGREE_Y_${DEGREE_Y}_${BSPLINES_TYPE}_${EVALUATOR}_${BCL}_${BCR}")
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
