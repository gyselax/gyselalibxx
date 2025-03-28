

add_subdirectory(connectivity)
add_subdirectory(data_types)
add_subdirectory(spline)
add_subdirectory(utils)

if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/interface_derivatives/CMakeLists.txt)
add_subdirectory(interface_derivatives)
endif()
