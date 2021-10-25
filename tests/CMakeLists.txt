cmake_minimum_required(VERSION 3.15)

include(GoogleTest)

add_executable(unit_tests
    fft.cpp
    splitvlasovsolver.cpp
)
target_compile_features(unit_tests PUBLIC cxx_std_17)
target_link_libraries(unit_tests
    PUBLIC
        GTest::gtest_main
        GTest::gmock
        GTest::gmock_main
        FFTW::Double
        paraconf::paraconf
        DDC::PDI
        vcx::fft
        vcx::vlasov
)
gtest_discover_tests(unit_tests)

add_subdirectory(landau)
