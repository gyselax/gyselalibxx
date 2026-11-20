

set(test_name_gtest "multipatch_interpolation_tests")
add_executable("${test_name_gtest}"
        interpolation_builder_2d_collection.cpp
        interpolation_builder_collection.cpp
        ../../main.cpp
)
target_link_libraries("${test_name_gtest}"
PUBLIC
        DDC::core
        GTest::gtest
        GTest::gmock

        gyselalibxx::multipatch_geometries
        gyselalibxx::multipatch_interpolation
        gyselalibxx::utils
)
target_compile_definitions("${test_name_gtest}" PUBLIC)
gtest_discover_tests("${test_name_gtest}" DISCOVERY_TIMEOUT ${GYSELALIBXX_DISCOVERY_TIMEOUT} DISCOVERY_MODE PRE_TEST)
