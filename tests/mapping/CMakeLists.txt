include(GoogleTest)

add_executable(unit_tests_mapping
    ../main.cpp
    mapping_analytical_inverse.cpp
    mapping_execution_space_access.cpp
    mapping_jacobian.cpp
    mapping_jacobian_matrix_coef.cpp
    mapping_static_properties.cpp
    metric_tensor_evaluator.cpp
    pseudo_cartesian_jacobian_matrix.cpp
    refined_discrete_mapping.cpp
)
target_compile_features(unit_tests_mapping PUBLIC cxx_std_17)
target_link_libraries(unit_tests_mapping
    PUBLIC
        GTest::gtest
        GTest::gmock
        gslx::mapping
        gslx::utils
)
gtest_discover_tests(unit_tests_mapping DISCOVERY_MODE PRE_TEST)
