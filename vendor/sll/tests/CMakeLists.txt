cmake_minimum_required(VERSION 3.15)

include(GoogleTest)

add_executable(splines_tests
    chunk.cpp
    periodic_spline_builder.cpp
    non_periodic_spline_builder.cpp
    view.cpp
    gauss_legendre_integration.cpp
)
if(VOICEXX_ENABLE_DEPRECATED)
    target_sources(splines_tests
	PUBLIC
        bsplines_non_uniform.cpp
        bsplines_uniform.cpp
        splines.cpp
    )
endif()
target_compile_features(splines_tests PUBLIC cxx_std_17)
target_link_libraries(splines_tests
    PUBLIC
        GTest::gtest_main
        sll::splines
)
gtest_discover_tests(splines_tests)

add_executable(bsplines_tests
	bsplines.cpp
)
target_compile_features(bsplines_tests PUBLIC cxx_std_17)
target_link_libraries(bsplines_tests
    PUBLIC
        GTest::gtest_main
        sll::splines
)

gtest_discover_tests(bsplines_tests)

add_executable(matrix_tests
	matrix.cpp
)
target_compile_features(matrix_tests PUBLIC cxx_std_17)
target_link_libraries(matrix_tests
    PUBLIC
        GTest::gtest_main
        sll::splines
)
gtest_discover_tests(matrix_tests)
