add_library(math_tools INTERFACE)

target_compile_features(math_tools
	INTERFACE
		cxx_std_17
)

target_include_directories(math_tools
    INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries(math_tools
	INTERFACE
		vcx::array_tools
)

add_library(sll::math_tools ALIAS math_tools)
