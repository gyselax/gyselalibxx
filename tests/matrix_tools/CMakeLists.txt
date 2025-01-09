add_executable(matrix_tests
  ../main.cpp
  matrix.cpp
  matrix_batch_ell.cpp
  matrix_batch_csr.cpp
  matrix_batch_tridiag.cpp
)
target_link_libraries(matrix_tests
    PUBLIC
        GTest::gtest
        GTest::gmock
        sll::SLL
        gslx::matrix_tools
)
gtest_discover_tests(matrix_tests DISCOVERY_MODE PRE_TEST)

