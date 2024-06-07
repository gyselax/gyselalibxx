# SPDX-License-Identifier: MIT

cmake_minimum_required(VERSION 3.15)

include(GoogleTest)

add_executable(unit_tests_pde_solvers
    fftpoissonsolver.cpp
    ../main.cpp
)
target_compile_features(unit_tests_pde_solvers PUBLIC cxx_std_17)
target_link_libraries(unit_tests_pde_solvers
    PUBLIC
        DDC::DDC
        GTest::gtest
        GTest::gmock
        gslx::utils
        gslx::pde_solvers
)

gtest_discover_tests(unit_tests_pde_solvers)
