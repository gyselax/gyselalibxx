add_library(array_tools INTERFACE)

target_include_directories(array_tools
	INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries(array_tools
	INTERFACE
		std::mdspan
		std::mdarray
)

add_library(vcx::array_tools ALIAS array_tools)
