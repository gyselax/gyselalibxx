include(GoogleTest)

add_executable(misc_tests_sll
    main.cpp
    view.cpp
    gauss_legendre_integration.cpp
)
target_compile_features(misc_tests_sll PUBLIC cxx_std_17)
target_link_libraries(misc_tests_sll
    PUBLIC
        GTest::gtest
        sll::SLL
)

add_executable(bernstein_tests_sll
  main.cpp
	bernstein.cpp
)
target_compile_features(bernstein_tests_sll PUBLIC cxx_std_17)
target_link_libraries(bernstein_tests_sll
    PUBLIC
        GTest::gtest
        sll::SLL
)
gtest_discover_tests(bernstein_tests_sll DISCOVERY_MODE PRE_TEST)

add_executable(metric_tensor_tests_sll
  main.cpp
	metric_tensor.cpp
)
target_compile_features(metric_tensor_tests_sll PUBLIC cxx_std_17)
target_link_libraries(metric_tensor_tests_sll
    PUBLIC
        GTest::gtest
        sll::SLL
)
gtest_discover_tests(metric_tensor_tests_sll DISCOVERY_MODE PRE_TEST)

add_executable(jacobian_mapping_tests_sll
  main.cpp
	mapping_jacobian.cpp
)
target_compile_features(jacobian_mapping_tests_sll PUBLIC cxx_std_17)
target_link_libraries(jacobian_mapping_tests_sll
    PUBLIC
        GTest::gtest
        sll::SLL
)
gtest_discover_tests(jacobian_mapping_tests_sll DISCOVERY_MODE PRE_TEST)

add_executable(jacobian_mapping_matrix_coef_tests_sll
  main.cpp
	mapping_jacobian_matrix_coef.cpp
)
target_compile_features(jacobian_mapping_matrix_coef_tests_sll PUBLIC cxx_std_17)
target_link_libraries(jacobian_mapping_matrix_coef_tests_sll
    PUBLIC
        GTest::gtest
        sll::SLL
)
gtest_discover_tests(jacobian_mapping_matrix_coef_tests_sll DISCOVERY_MODE PRE_TEST)

add_executable(pseudo_cartesian_tests_sll
  main.cpp
	pseudo_cartesian_jacobian_matrix.cpp
)
target_compile_features(pseudo_cartesian_tests_sll PUBLIC cxx_std_17)
target_link_libraries(pseudo_cartesian_tests_sll
    PUBLIC
        GTest::gtest
        sll::SLL
)
gtest_discover_tests(pseudo_cartesian_tests_sll DISCOVERY_MODE PRE_TEST)

add_executable(refined_discrete_mapping_test
  main.cpp
	refined_discrete_mapping.cpp
)
target_compile_features(refined_discrete_mapping_test PUBLIC cxx_std_17)
target_link_libraries(refined_discrete_mapping_test
    PUBLIC
        GTest::gtest
        sll::SLL
)
gtest_discover_tests(refined_discrete_mapping_test DISCOVERY_MODE PRE_TEST)

add_executable(mapping_execution_space_access
  main.cpp
	mapping_execution_space_access.cpp
)
target_compile_features(mapping_execution_space_access PUBLIC cxx_std_17)
target_link_libraries(mapping_execution_space_access
    PUBLIC
        GTest::gtest
        sll::SLL
)
gtest_discover_tests(mapping_execution_space_access DISCOVERY_MODE PRE_TEST)

