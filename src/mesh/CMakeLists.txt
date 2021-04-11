add_library(mesh INTERFACE)

target_compile_features(mesh
	INTERFACE
		cxx_std_17
)

target_include_directories(mesh
	INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries(mesh
	INTERFACE
		vcx::array_tools
		vcx::geometry
)

add_library(vcx::mesh ALIAS mesh)
