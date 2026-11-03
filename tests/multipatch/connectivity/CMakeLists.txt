

add_executable("unit_tests_connectivity"
        coord_transformation.cpp
        coord_transformation_periodic.cpp
        index_transformation_uniform.cpp
        index_transformation_non_uniform.cpp
        index_transformation_uniform_periodic.cpp
        index_transformation_non_uniform_periodic.cpp
        matching_idx_slice_uniform_grids.cpp
        matching_idx_slice_non_uniform_grids.cpp
        multipatch_connectivity.cpp
        patch_locator_onion_shape_2patches.cpp
        ../../main.cpp
)
target_link_libraries("unit_tests_connectivity"
PUBLIC
        DDC::core
        GTest::gtest
        GTest::gmock

        gyselalibxx::multipatch_connectivity
        gyselalibxx::multipatch_data_types
        gyselalibxx::multipatch_geometries
        gyselalibxx::utils
)
target_compile_definitions("unit_tests_connectivity" PUBLIC)
gtest_discover_tests("unit_tests_connectivity" DISCOVERY_TIMEOUT ${GYSELALIBXX_DISCOVERY_TIMEOUT} DISCOVERY_MODE PRE_TEST)
