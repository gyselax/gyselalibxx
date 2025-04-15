set(test_name_gtest "Ginkgo_test")
add_executable("${test_name_gtest}"
    test_with_ginkgo.cpp
    #../../main.cpp
)
target_link_libraries("${test_name_gtest}"
PUBLIC
    DDC::core
    gslx::utils
    # GTest::gtest
    # GTest::gmock
)
target_compile_definitions("${test_name_gtest}" PUBLIC)
gtest_discover_tests("${test_name_gtest}" DISCOVERY_MODE PRE_TEST)
