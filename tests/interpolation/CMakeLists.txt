# SPDX-License-Identifier: MIT

if (DEFINED ENV{POLAR_SPLINES_TEST_DEGREE_MIN})
    set(POLAR_SPLINES_TEST_DEGREE_MIN $ENV{POLAR_SPLINES_TEST_DEGREE_MIN})
else()
    set(POLAR_SPLINES_TEST_DEGREE_MIN 3 CACHE STRING "Minimum degree to test splines.")
endif()

if (DEFINED ENV{POLAR_SPLINES_TEST_DEGREE_MAX})
    set(POLAR_SPLINES_TEST_DEGREE_MAX $ENV{POLAR_SPLINES_TEST_DEGREE_MAX})
else()
    set(POLAR_SPLINES_TEST_DEGREE_MAX 3 CACHE STRING "Maximum degree to test splines.")
endif()

include(GoogleTest)

add_executable(bernstein_tests
    ../main.cpp
    bernstein.cpp
)
target_compile_features(bernstein_tests PUBLIC cxx_std_17)
target_link_libraries(bernstein_tests
    PUBLIC
        GTest::gtest
        GTest::gmock
        gslx::polar_splines
)
target_include_directories(bernstein_tests
    PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}/.."
)
gtest_discover_tests(bernstein_tests DISCOVERY_MODE PRE_TEST)

add_executable(polar_bsplines_tests
    ../main.cpp
    polar_bsplines.cpp
    polar_spline_create_mirror.cpp
)
target_compile_features(polar_bsplines_tests PUBLIC cxx_std_17)
target_link_libraries(polar_bsplines_tests
    PUBLIC
        GTest::gtest
        GTest::gmock
        gslx::polar_splines
)
target_include_directories(polar_bsplines_tests
    PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}/.."
)

gtest_discover_tests(polar_bsplines_tests DISCOVERY_MODE PRE_TEST)

foreach(CONTINUITY RANGE -1 1)
  math(EXPR MIN_DEGREE "${CONTINUITY}+1")
  if (${MIN_DEGREE} LESS 1)
    set(MIN_DEGREE 1)
  endif()
  if (${MIN_DEGREE} LESS ${POLAR_SPLINES_TEST_DEGREE_MIN})
    set(MIN_DEGREE ${POLAR_SPLINES_TEST_DEGREE_MIN})
  endif()
  foreach(DEGREE_R RANGE ${MIN_DEGREE} ${POLAR_SPLINES_TEST_DEGREE_MAX})
    foreach(DEGREE_P RANGE ${MIN_DEGREE} ${POLAR_SPLINES_TEST_DEGREE_MAX})
      foreach(BSPLINES_TYPE "BSPLINES_TYPE_UNIFORM" "BSPLINES_TYPE_NON_UNIFORM")
        foreach(MAPPING_TYPE "CIRCULAR_MAPPING" "CZARNY_MAPPING")
          set(test_name "polar_splines_tests_DEGREE_R_${DEGREE_R}_DEGREE_P_${DEGREE_P}_CONTINUITY_${CONTINUITY}_${BSPLINES_TYPE}_${MAPPING_TYPE}")
          add_executable("${test_name}" polar_splines.cpp)
          target_compile_features("${test_name}" PUBLIC cxx_std_17)
          target_link_libraries("${test_name}"
            PUBLIC
            GTest::gtest
            GTest::gmock
            gslx::polar_splines
          )
          target_compile_definitions("${test_name}" PUBLIC -DDEGREE_R=${DEGREE_R} -DDEGREE_P=${DEGREE_P} -D${BSPLINES_TYPE} -DCONTINUITY=${CONTINUITY} -D${MAPPING_TYPE})
          gtest_discover_tests("${test_name}"
              DISCOVERY_MODE PRE_TEST
              TEST_SUFFIX "_${test_name}")
        endforeach()
      endforeach()
    endforeach()
  endforeach()
endforeach()

