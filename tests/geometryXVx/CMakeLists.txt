# SPDX-License-Identifier: MIT

cmake_minimum_required(VERSION 3.15)

include(GoogleTest)

foreach(GEOMETRY_VARIANT IN LISTS GEOMETRY_XVx_VARIANTS_LIST)

add_executable(unit_tests_${GEOMETRY_VARIANT}
    collisions_intra_gridvx.cpp
    collisions_intra_maxwellian.cpp
    fluid_moments.cpp
    kineticsource.cpp
    krooksource.cpp
    main.cpp
    masks.cpp
    quadrature.cpp
    species_info.cpp
    splitvlasovsolver.cpp
    rk2.cpp
)
target_compile_features(unit_tests_${GEOMETRY_VARIANT} PUBLIC cxx_std_17)
target_link_libraries(unit_tests_${GEOMETRY_VARIANT}
    PUBLIC
        DDC::PDI_Wrapper
        GTest::gtest
        GTest::gmock
        GTest::gmock_main
        paraconf::paraconf
        sll::splines
        vcx::advection
        vcx::boltzmann_${GEOMETRY_VARIANT}
        vcx::initialization_${GEOMETRY_VARIANT}
        vcx::poisson_${GEOMETRY_VARIANT}
        vcx::quadrature
        vcx::utils_${GEOMETRY_VARIANT}
)

gtest_discover_tests(unit_tests_${GEOMETRY_VARIANT}
    TEST_SUFFIX "_${GEOMETRY_VARIANT}")

endforeach()

target_sources(unit_tests_xperiod_vx
    PRIVATE
        fft.cpp
        femperiodicpoissonsolver.cpp
)

target_link_libraries(unit_tests_xperiod_vx PRIVATE vcx::fft)

target_sources(unit_tests_xnonperiod_vx PRIVATE femnonperiodicpoissonsolver.cpp)

add_subdirectory(bump_on_tail)

add_subdirectory(landau)
