




set(test_name_gtest "interface")

add_executable("${test_name_gtest}"
        coord_transformation.cpp
        ../main.cpp
)
target_link_libraries("${test_name_gtest}"
PUBLIC
        DDC::DDC
        GTest::gtest
        GTest::gmock

        gslx::multipatch_interface
)
target_compile_definitions("${test_name_gtest}" PUBLIC)
gtest_discover_tests("${test_name_gtest}" )
