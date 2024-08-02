

add_executable("coord_transformation_test"
        coord_transformation.cpp
        ../../main.cpp
)
target_link_libraries("coord_transformation_test"
PUBLIC
        DDC::DDC
        GTest::gtest
        GTest::gmock

        gslx::multipatch_connectivity
        gslx::multipatch_geometries
)
target_compile_definitions("coord_transformation_test" PUBLIC)
gtest_discover_tests("coord_transformation_test" DISCOVERY_MODE PRE_TEST)



add_executable("index_transformation_test"
        index_transformation_uniform.cpp
        index_transformation_non_uniform.cpp
        ../../main.cpp
)
target_link_libraries("index_transformation_test"
PUBLIC
        DDC::DDC
        GTest::gtest
        GTest::gmock

        gslx::multipatch_connectivity
        gslx::multipatch_geometries
)
target_compile_definitions("index_transformation_test" PUBLIC)
gtest_discover_tests("index_transformation_test" DISCOVERY_MODE PRE_TEST)

