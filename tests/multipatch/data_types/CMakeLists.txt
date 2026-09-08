
set(test_name_gtest "multipatch_field_test")
add_executable("${test_name_gtest}"
        multipatch_field_2p.cpp
        multipatch_field_9p.cpp
        ../../main.cpp
)
target_link_libraries("${test_name_gtest}"
PUBLIC
        DDC::core
        GTest::gtest
        GTest::gmock

        gyselalibxx::multipatch_data_types
        gyselalibxx::utils
        gyselalibxx::multipatch_connectivity
        gyselalibxx::multipatch_geometries
)
target_compile_definitions("${test_name_gtest}" PUBLIC)
gtest_discover_tests("${test_name_gtest}" DISCOVERY_TIMEOUT ${GYSELALIBXX_DISCOVERY_TIMEOUT} DISCOVERY_MODE PRE_TEST)
