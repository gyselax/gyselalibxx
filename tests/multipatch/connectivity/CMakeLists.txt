

add_executable("unit_tests_connectivity"
        coord_transformation.cpp
        index_transformation_uniform.cpp
        index_transformation_non_uniform.cpp
        multipatch_connectivity.cpp
        matching_idx_slice_uniform_grids.cpp
        matching_idx_slice_non_uniform_grids.cpp
        ../../main.cpp
)
target_link_libraries("unit_tests_connectivity"
PUBLIC
        DDC::DDC
        GTest::gtest
        GTest::gmock

        gslx::multipatch_connectivity
        gslx::multipatch_geometries
        gslx::utils
)
target_compile_definitions("unit_tests_connectivity" PUBLIC)
gtest_discover_tests("unit_tests_connectivity" DISCOVERY_MODE PRE_TEST)
