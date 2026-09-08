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
        DDC::core
        DDC::splines
        gyselalibxx::matrix_tools
        gyselalibxx::math_tools
)
gtest_discover_tests(matrix_tests DISCOVERY_TIMEOUT ${GYSELALIBXX_DISCOVERY_TIMEOUT} DISCOVERY_MODE PRE_TEST)

