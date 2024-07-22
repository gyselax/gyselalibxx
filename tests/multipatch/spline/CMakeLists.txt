

set(test_name_gtest "multipatch_spline_builder")
add_executable("${test_name_gtest}"
        multipatch_spline_builder.cpp
        ../../main.cpp
)
target_link_libraries("${test_name_gtest}"
PUBLIC
        DDC::DDC
        GTest::gtest
        GTest::gmock

        gslx::multipatch_spline
        gslx::utils
        
)
target_compile_definitions("${test_name_gtest}" PUBLIC)
gtest_discover_tests("${test_name_gtest}" DISCOVERY_MODE PRE_TEST)
