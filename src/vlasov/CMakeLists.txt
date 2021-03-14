add_library(vlasov STATIC
	advection1d.cpp
	vlasov.cpp
)

target_include_directories( vlasov
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries( vlasov
	PUBLIC
		vcx::array_tools
		sll::splines
		vcx::mesh
)

add_library( vcx::vlasov ALIAS vlasov )
