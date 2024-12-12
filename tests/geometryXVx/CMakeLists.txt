# SPDX-License-Identifier: MIT

include(GoogleTest)

set(BASIC_GEOMETRY_XVx_VARIANTS_LIST xperiod_vx xnonperiod_vx)

foreach(GEOMETRY_VARIANT IN LISTS BASIC_GEOMETRY_XVx_VARIANTS_LIST)

add_executable(unit_tests_${GEOMETRY_VARIANT}
    collisions_inter.cpp
    collisions_intra_gridvx.cpp
    collisions_intra_maxwellian.cpp
    fluid_moments.cpp
    kineticsource.cpp
    krooksource.cpp
    masks.cpp
    splitvlasovsolver.cpp
    maxwellian.cpp
    ../main.cpp
)

target_link_libraries(unit_tests_${GEOMETRY_VARIANT}
    PUBLIC
        DDC::pdi
        GTest::gtest
        GTest::gmock
        paraconf::paraconf
        gslx::advection
        gslx::boltzmann_${GEOMETRY_VARIANT}
        gslx::moments
        gslx::poisson_${GEOMETRY_VARIANT}
        gslx::quadrature
        gslx::speciesinfo
        gslx::time_integration_${GEOMETRY_VARIANT}
        gslx::utils_${GEOMETRY_VARIANT}
)

gtest_discover_tests(unit_tests_${GEOMETRY_VARIANT}
    TEST_SUFFIX "_${GEOMETRY_VARIANT}"
    PROPERTIES TIMEOUT 10
    DISCOVERY_MODE PRE_TEST
)

endforeach()

add_executable(unit_tests_lagrange
lagrange_test.cpp
    ../main.cpp
)

target_link_libraries(unit_tests_lagrange
    PUBLIC
        GTest::gtest
        GTest::gmock
        gslx::interpolation
        gslx::utils


)
gtest_discover_tests(unit_tests_lagrange
    TEST_SUFFIX "_${lagrange}"
    DISCOVERY_MODE PRE_TEST)

add_executable(unit_tests_advection
    velocityadvection.cpp
    spatialadvection.cpp
    ../main.cpp
)

target_link_libraries(unit_tests_advection	 
    PUBLIC
        GTest::gtest
        GTest::gmock
        gslx::interpolation
	gslx::advection
	gslx::geometry_xperiod_vx
	gslx::utils


)
gtest_discover_tests(unit_tests_advection 
    TEST_SUFFIX "_${advection}"
    DISCOVERY_MODE PRE_TEST)


# Set variables for restart tests
if (DEFINED ENV{RELATIVE_RESTART_TOLERANCE})
    set(RELATIVE_RESTART_TOLERANCE $ENV{RELATIVE_RESTART_TOLERANCE})
else()
    set(RELATIVE_RESTART_TOLERANCE "1e-14")
endif()
if (DEFINED ENV{ABSOLUTE_RESTART_TOLERANCE})
    set(ABSOLUTE_RESTART_TOLERANCE $ENV{ABSOLUTE_RESTART_TOLERANCE})
else()
    set(ABSOLUTE_RESTART_TOLERANCE "1e-14")
endif()

add_subdirectory(bump_on_tail)

add_subdirectory(landau)
