# SPDX-License-Identifier: MIT


include(GoogleTest)

foreach(GEOMETRY_VARIANT IN LISTS GEOMETRY_XVx_VARIANTS_LIST)

add_executable(unit_tests_geometryMX_${GEOMETRY_VARIANT}
	../../main.cpp
    fluidspecies.cpp
    moments.cpp
    predcorr_hybrid.cpp
    diffusiveneutrals.cpp
    ratecoefficients.cpp
)

target_link_libraries(unit_tests_geometryMX_${GEOMETRY_VARIANT}
    PUBLIC
        GTest::gtest
        GTest::gmock
        gslx::advection
        gslx::boltzmann_${GEOMETRY_VARIANT}
        gslx::fluidinitialization_${GEOMETRY_VARIANT}
        gslx::fluidsolver_${GEOMETRY_VARIANT}
        gslx::moments
        gslx::poisson_${GEOMETRY_VARIANT}
        gslx::quadrature
        gslx::speciesinfo
        gslx::time_integration_${GEOMETRY_VARIANT}
        gslx::time_integration_hybrid_${GEOMETRY_VARIANT}
        gslx::utils_${GEOMETRY_VARIANT}
)

gtest_discover_tests(unit_tests_geometryMX_${GEOMETRY_VARIANT}
    TEST_SUFFIX "_${GEOMETRY_VARIANT}"
    PROPERTIES TIMEOUT 10
    DISCOVERY_MODE PRE_TEST
)

endforeach()
