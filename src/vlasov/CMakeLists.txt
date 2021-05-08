add_library(vlasov STATIC
	splitvlasovsolver.cpp
	vlasovsolver.cpp
	advection1d.cpp
)

target_compile_features(vlasov
	PUBLIC
		cxx_std_17
)

target_include_directories(vlasov
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries(vlasov
	PUBLIC
		vcx::array_tools
		sll::splines
		vcx::mesh
)

add_library(vcx::vlasov ALIAS vlasov)
