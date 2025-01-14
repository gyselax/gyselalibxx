include(GoogleTest)

add_executable(misc_tests_sll
    main.cpp
    view.cpp
)
target_compile_features(misc_tests_sll PUBLIC cxx_std_17)
target_link_libraries(misc_tests_sll
    PUBLIC
        GTest::gtest
        sll::SLL
)
