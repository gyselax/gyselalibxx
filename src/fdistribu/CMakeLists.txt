cmake_minimum_required(VERSION 3.15)

add_library(fdistribu fdistribu.cpp)
target_include_directories(fdistribu
	PUBLIC
		"$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
)
target_link_libraries(fdistribu PUBLIC
	DDC::DDC
	sll::splines
	vcx::geometry
)
add_library(vcx::fdistribu ALIAS fdistribu)