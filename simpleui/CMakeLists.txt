cmake_minimum_required(VERSION 3.18)

add_library(simpleui
		src/commandlineconfig.cpp
		src/screenprinter.cpp
)
target_include_directories(simpleui PUBLIC include/)
target_link_libraries(simpleui PUBLIC GLCS2020_project::baselib)
add_library(GLCS2020_project::simpleui ALIAS simpleui)
