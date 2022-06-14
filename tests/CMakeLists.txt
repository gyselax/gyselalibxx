# SPDX-License-Identifier: MIT

cmake_minimum_required(VERSION 3.15)

include(GoogleTest)

add_executable(unit_tests main.cpp splitvlasovsolver.cpp)
target_compile_features(unit_tests PUBLIC cxx_std_17)
target_link_libraries(unit_tests
    PUBLIC
        GTest::gtest
        GTest::gmock
        GTest::gmock_main
        paraconf::paraconf
        DDC::PDI_Wrapper
        vcx::vlasov
        vcx::poisson
        vcx::quadrature
        sll::splines
)

add_subdirectory(bump_on_tail)
if(VOICEXX_ENABLE_PERIODIC_RDIMX)
    target_sources(unit_tests
        PRIVATE
            fft.cpp
            femperiodicpoissonsolver.cpp
    )

    target_link_libraries(unit_tests PRIVATE vcx::fft)

    add_subdirectory(landau)
else()
    target_sources(unit_tests PRIVATE femnonperiodicpoissonsolver.cpp)
endif()
target_sources(unit_tests PRIVATE quadrature.cpp)

gtest_discover_tests(unit_tests)
