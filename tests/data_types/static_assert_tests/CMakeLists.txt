
function(add_static_assert_test file error_msg)

add_executable(${file}-build-fail-test
    ${file}.cpp
)
target_link_libraries(${file}-build-fail-test
        PUBLIC
        gslx::data_types
        gslx::utils
)

set_target_properties(${file}-build-fail-test PROPERTIES
    EXCLUDE_FROM_ALL TRUE
    EXCLUDE_FROM_DEFAULT_BUILD TRUE
)

add_test(
    NAME ${file}
    COMMAND ${CMAKE_COMMAND} --build . --target ${file}-build-fail-test --config $<CONFIG>
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR} # the place where we would run ctest
)

set_tests_properties(${file} PROPERTIES
    PASS_REGULAR_EXPRESSION ${error_msg}
)

set_property(TEST ${file} PROPERTY TIMEOUT 20)

endfunction()

add_static_assert_test(tensor_mul_cov_cov "Repeated indices should not be associated with two covariant or two contravariant indices.")
add_static_assert_test(tensor_mul_contra_contra "Repeated indices should not be associated with two covariant or two contravariant indices.")
add_static_assert_test(tensor_mul_incompatible_index_set "Cannot sum over incompatible VectorIndexSets.")
add_static_assert_test(tensor_mul_3_indices "You should not have more than two of any one index in an index expression.")
add_static_assert_test(tensor_mul_multiple_types "All tensors in a tensor_mul must have the same element type.")
add_static_assert_test(tensor_mul_wrong_index "Index is not compatible with tensor type")


