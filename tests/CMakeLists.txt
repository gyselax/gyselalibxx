# SPDX-License-Identifier: MIT

cmake_minimum_required(VERSION 3.15)

include(GoogleTest)

foreach(GEOMETRY_VARIANT IN LISTS GEOMETRY_VARIANTS_LIST)

add_executable(unit_tests_${GEOMETRY_VARIANT}
    main.cpp
    quadrature.cpp
    splitvlasovsolver.cpp
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
        vcx::poisson_${GEOMETRY_VARIANT}
        vcx::quadrature
        vcx::vlasov_${GEOMETRY_VARIANT}
)

gtest_discover_tests(unit_tests_${GEOMETRY_VARIANT}
    TEST_SUFFIX "_${GEOMETRY_VARIANT}")

endforeach()

target_sources(unit_tests_periodic_rdimx
    PRIVATE
        fft.cpp
        femperiodicpoissonsolver.cpp
)

target_link_libraries(unit_tests_periodic_rdimx PRIVATE vcx::fft)

target_sources(unit_tests_nonperiodic_rdimx PRIVATE femnonperiodicpoissonsolver.cpp)

add_subdirectory(bump_on_tail)

add_subdirectory(landau)
