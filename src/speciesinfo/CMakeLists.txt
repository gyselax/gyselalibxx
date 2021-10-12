cmake_minimum_required(VERSION 3.15)

add_library(speciesinfo species_info.cpp)
target_include_directories(speciesinfo
	PUBLIC
		"$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
)
target_link_libraries(speciesinfo PUBLIC
	DDC::DDC
	sll::splines
	vcx::geometry
)
add_library(vcx::speciesinfo ALIAS speciesinfo)
