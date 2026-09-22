# SPDX-License-Identifier: MIT

if (DEFINED ENV{GYSELALIBXX_POLAR_SPLINES_TEST_DEGREE_MIN})
    set(GYSELALIBXX_POLAR_SPLINES_TEST_DEGREE_MIN $ENV{GYSELALIBXX_POLAR_SPLINES_TEST_DEGREE_MIN})
else()
    set(GYSELALIBXX_POLAR_SPLINES_TEST_DEGREE_MIN 3 CACHE STRING "Minimum degree to test splines.")
endif()

if (DEFINED ENV{GYSELALIBXX_POLAR_SPLINES_TEST_DEGREE_MAX})
    set(GYSELALIBXX_POLAR_SPLINES_TEST_DEGREE_MAX $ENV{GYSELALIBXX_POLAR_SPLINES_TEST_DEGREE_MAX})
else()
    set(GYSELALIBXX_POLAR_SPLINES_TEST_DEGREE_MAX 3 CACHE STRING "Maximum degree to test splines.")
endif()

include(GoogleTest)

add_executable(bernstein_tests
    ../main.cpp
    bernstein.cpp
)
target_link_libraries(bernstein_tests
    PUBLIC
        GTest::gtest
        GTest::gmock
        gyselalibxx::polar_splines
)
target_include_directories(bernstein_tests
    PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}/.."
)
gtest_discover_tests(bernstein_tests DISCOVERY_TIMEOUT ${GYSELALIBXX_DISCOVERY_TIMEOUT} DISCOVERY_MODE PRE_TEST)

add_executable(polar_bsplines_tests
    ../main.cpp
    polar_bsplines.cpp
)
target_link_libraries(polar_bsplines_tests
    PUBLIC
        GTest::gtest
        GTest::gmock
        gyselalibxx::polar_splines
)
target_include_directories(polar_bsplines_tests
    PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}/.."
)

gtest_discover_tests(polar_bsplines_tests DISCOVERY_TIMEOUT ${GYSELALIBXX_DISCOVERY_TIMEOUT} DISCOVERY_MODE PRE_TEST)

foreach(CONTINUITY RANGE -1 1)
  math(EXPR MIN_DEGREE "${CONTINUITY}+1")
  if (${MIN_DEGREE} LESS 1)
    set(MIN_DEGREE 1)
  endif()
  if (${MIN_DEGREE} LESS ${GYSELALIBXX_POLAR_SPLINES_TEST_DEGREE_MIN})
    set(MIN_DEGREE ${GYSELALIBXX_POLAR_SPLINES_TEST_DEGREE_MIN})
  endif()
  foreach(DEGREE_R RANGE ${MIN_DEGREE} ${GYSELALIBXX_POLAR_SPLINES_TEST_DEGREE_MAX})
    foreach(DEGREE_P RANGE ${MIN_DEGREE} ${GYSELALIBXX_POLAR_SPLINES_TEST_DEGREE_MAX})
      foreach(BSPLINES_TYPE "BSPLINES_TYPE_UNIFORM" "BSPLINES_TYPE_NON_UNIFORM")
        foreach(MAPPING_TYPE "CIRCULAR_MAPPING" "CZARNY_MAPPING")
          set(test_name "polar_splines_tests_DEGREE_R_${DEGREE_R}_DEGREE_P_${DEGREE_P}_CONTINUITY_${CONTINUITY}_${BSPLINES_TYPE}_${MAPPING_TYPE}")
          add_executable("${test_name}" polar_splines.cpp)
          target_link_libraries("${test_name}"
            PUBLIC
            GTest::gtest
            GTest::gmock
            gyselalibxx::polar_splines
          )
          target_compile_definitions("${test_name}" PUBLIC -DDEGREE_R=${DEGREE_R} -DDEGREE_P=${DEGREE_P} -D${BSPLINES_TYPE} -DCONTINUITY=${CONTINUITY} -D${MAPPING_TYPE})
          gtest_discover_tests("${test_name}"
              TEST_SUFFIX "_${test_name}"
              DISCOVERY_TIMEOUT ${GYSELALIBXX_DISCOVERY_TIMEOUT}
              DISCOVERY_MODE PRE_TEST)
        endforeach()
      endforeach()
    endforeach()
  endforeach()
endforeach()


add_executable(lagrange_tests
    ../main.cpp
    lagrange_bsplines.cpp
    lagrange_interpolation.cpp
    nd_lagrange_interpolation.cpp
)
target_link_libraries(lagrange_tests
    PUBLIC
        GTest::gtest
        GTest::gmock
        gyselalibxx::interpolation
        gyselalibxx::math_tools
)
target_include_directories(lagrange_tests
    PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}/.."
)

gtest_discover_tests(lagrange_tests DISCOVERY_TIMEOUT ${GYSELALIBXX_DISCOVERY_TIMEOUT} DISCOVERY_MODE PRE_TEST)

add_executable(spline_interpolation_tests
    ../main.cpp
    spline_interpolation.cpp
)
target_link_libraries(spline_interpolation_tests
    PUBLIC
        GTest::gtest
        GTest::gmock
        gyselalibxx::interpolation
)

gtest_discover_tests(spline_interpolation_tests DISCOVERY_TIMEOUT ${GYSELALIBXX_DISCOVERY_TIMEOUT} DISCOVERY_MODE PRE_TEST)

# MPISplineBuilderLocalCrouseilles needs multiple MPI ranks, so (unlike the executables above)
# it is not registered with gtest_discover_tests but with individual mpirun-wrapped ctest
# entries, following the pattern in tests/mpi_parallelisation/CMakeLists.txt.
add_executable(mpi_spline_builder_local_crouseilles_tests
    mpi_spline_builder_local_crouseilles.cpp
)
target_link_libraries(mpi_spline_builder_local_crouseilles_tests
    PUBLIC
        DDC::core
        GTest::gtest
        GTest::gmock
        gyselalibxx::interpolation
        gyselalibxx::mpi_parallelisation
)

function(make_mpi_spline_builder_local_crouseilles_test test_name)
    add_test(NAME ${test_name}
        COMMAND
        "${MPIEXEC_EXECUTABLE}"
        "-n"
        "2"
        "$<TARGET_FILE:mpi_spline_builder_local_crouseilles_tests>"
        "--gtest_filter=${test_name}"
    )
endfunction()

make_mpi_spline_builder_local_crouseilles_test(MPISplineBuilderLocalCrouseilles.ReproducesLinearFunction)
make_mpi_spline_builder_local_crouseilles_test(MPISplineBuilderLocalCrouseilles.ConvergesForCosine)

