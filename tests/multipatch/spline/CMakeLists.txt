

set(test_name_gtest "multipatch_spline_tests")
add_executable("${test_name_gtest}"
        constant_extrapolation_rules_onion.cpp
        multipatch_spline_evaluator.cpp
        spline_builder_collection.cpp
        spline_builder_2d_collection.cpp
        null_extrapolation_rules.cpp
        spline_testing_tools.cpp
        ../../main.cpp
)
target_link_libraries("${test_name_gtest}"
PUBLIC
        DDC::core
        GTest::gtest
        GTest::gmock

        gyselalibxx::multipatch_geometries
        gyselalibxx::multipatch_spline
        gyselalibxx::utils
)
target_compile_definitions("${test_name_gtest}" PUBLIC)
gtest_discover_tests("${test_name_gtest}" DISCOVERY_TIMEOUT ${GYSELALIBXX_DISCOVERY_TIMEOUT} DISCOVERY_MODE PRE_TEST)
