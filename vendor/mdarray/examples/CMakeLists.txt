
function(mdarray_add_example EXENAME)
  add_executable(${EXENAME} ${EXENAME}.cpp)
  target_link_libraries(${EXENAME} mdarray)
endfunction(mdarray_add_example)

add_subdirectory(dot_product)
