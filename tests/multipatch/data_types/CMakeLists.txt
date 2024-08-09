
set(test_name_gtest "multipatch_field_test")
add_executable("${test_name_gtest}"
        multipatch_field_2p.cpp
        multipatch_field_9p.cpp
        ../../main.cpp
)
target_link_libraries("${test_name_gtest}"
PUBLIC
        DDC::DDC
        GTest::gtest
        GTest::gmock

        gslx::multipatch_data_types
        gslx::utils
        gslx::multipatch_connectivity
        gslx::multipatch_geometries
)
target_compile_definitions("${test_name_gtest}" PUBLIC)
gtest_discover_tests("${test_name_gtest}" DISCOVERY_MODE PRE_TEST)
