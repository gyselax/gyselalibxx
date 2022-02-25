
macro(add_compilation_test name)
  add_executable(${name} ${name}.cpp)
  target_link_libraries(${name} mdarray)
endmacro()

add_compilation_test(ctest_constructor_sfinae)

