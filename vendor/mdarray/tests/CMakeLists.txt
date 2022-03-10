
macro(mdarray_add_test name)
  add_executable(${name} ${name}.cpp)
  target_link_libraries(${name} mdarray GTest::GTest GTest::Main)
  add_test(${name} ${name})
endmacro()

find_package(GTest REQUIRED)

#mdarray_add_test(test_thingy)

