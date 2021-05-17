cmake_minimum_required(VERSION 3.15)

add_library(mdsw INTERFACE)

target_include_directories(mdsw
	INTERFACE
		"$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
)

target_link_libraries(mdsw
	INTERFACE
		std::mdspan
)

add_library(mdspan_wrapper::mdsw ALIAS mdsw)

