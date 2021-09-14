cmake_minimum_required(VERSION 3.15)

add_library(geometry INTERFACE)
target_include_directories(geometry
	INTERFACE
		"$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
)
target_link_libraries(geometry INTERFACE
	ddc::ddc
	sll::splines
)
add_library(vcx::geometry ALIAS geometry)
