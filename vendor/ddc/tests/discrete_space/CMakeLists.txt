add_library(discrete_space_tests_lib STATIC discrete_space.cpp)
target_link_libraries(discrete_space_tests_lib
    PUBLIC
    GTest::gtest
    DDC::DDC
)

add_executable(discrete_space_tests main.cpp)
target_link_libraries(discrete_space_tests
    PUBLIC
        discrete_space_tests_lib
        GTest::gtest
        DDC::DDC
)

gtest_discover_tests(discrete_space_tests)
