

set(test_name_gtest "multipatch_interpolation_tests")
add_executable("${test_name_gtest}"
        interpolation_builder_2d_collection.cpp
        interpolation_builder_collection.cpp
        interpolation_evaluator_2d_collection.cpp
        ../spline/spline_testing_tools.cpp
        ../../main.cpp
)
target_include_directories("${test_name_gtest}"
PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/../spline
)
target_link_libraries("${test_name_gtest}"
PUBLIC
        DDC::core
        GTest::gtest
        GTest::gmock

        gyselalibxx::multipatch_geometries
        gyselalibxx::multipatch_interpolation
        gyselalibxx::multipatch_spline
        gyselalibxx::utils
)
target_compile_definitions("${test_name_gtest}" PUBLIC)
gtest_discover_tests("${test_name_gtest}" DISCOVERY_TIMEOUT ${GYSELALIBXX_DISCOVERY_TIMEOUT} DISCOVERY_MODE PRE_TEST)
