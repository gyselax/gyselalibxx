# SPDX-License-Identifier: MIT


include(GoogleTest)

add_executable(unit_tests_pde_solvers
    fftpoissonsolver.cpp
    femperiodicpoissonsolver.cpp
    femnonperiodicpoissonsolver.cpp
    ../main.cpp
)
target_link_libraries(unit_tests_pde_solvers
    PUBLIC
        DDC::core
        GTest::gtest
        GTest::gmock
        gslx::pde_solvers
        gslx::quadrature
        gslx::speciesinfo
        gslx::utils
)

gtest_discover_tests(unit_tests_pde_solvers DISCOVERY_MODE PRE_TEST)
