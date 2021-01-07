cmake_minimum_required(VERSION 3.18)

add_library(heatlib
		src/finitediffheatsolver.cpp
		src/fixedconditions.cpp
)
target_include_directories(heatlib PUBLIC include/)
target_link_libraries(heatlib PUBLIC GLCS2020_project::baselib)
add_library(GLCS2020_project::heatlib ALIAS heatlib)
