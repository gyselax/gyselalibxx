# SPDX-License-Identifier: MIT


include(GoogleTest)

add_executable(unit_tests_pde_solvers
    fftpoissonsolver.cpp
    ../main.cpp
)
target_link_libraries(unit_tests_pde_solvers
    PUBLIC
        DDC::DDC
        GTest::gtest
        GTest::gmock
        gslx::utils
        gslx::pde_solvers
)

gtest_discover_tests(unit_tests_pde_solvers)
