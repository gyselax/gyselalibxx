include(GoogleTest)

add_executable(unit_tests_coord_transformations
    ../main.cpp
    analytical_inverse.cpp
    coord_transformations_execution_space_access.cpp
    jacobian.cpp
    jacobian_matrix_coef.cpp
    coord_transformations_static_properties.cpp
    metric_tensor_evaluator.cpp
    pseudo_cartesian_jacobian_matrix.cpp
    refined_discrete_mapping.cpp
    vector_basis_change.cpp
)
target_link_libraries(unit_tests_coord_transformations
    PUBLIC
        GTest::gtest
        GTest::gmock
        gslx::coord_transformations
        gslx::speciesinfo
        gslx::utils
)
gtest_discover_tests(unit_tests_coord_transformations DISCOVERY_MODE PRE_TEST)
