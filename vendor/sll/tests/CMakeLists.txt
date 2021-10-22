cmake_minimum_required(VERSION 3.15)

include(GoogleTest)

add_executable(splines_tests
    chunk.cpp
    bsplines_non_uniform.cpp
    bsplines_uniform.cpp
    spline_builder.cpp
    splines.cpp
    view.cpp
    gauss_legendre_integration.cpp
)
target_compile_features(splines_tests PUBLIC cxx_std_17)
target_link_libraries(splines_tests
    PUBLIC
        GTest::gtest_main
        sll::splines
)
gtest_discover_tests(splines_tests)
