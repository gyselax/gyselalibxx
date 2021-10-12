add_library(vlasov STATIC
	nulladvectionvx.cpp
	splineadvectionx.cpp
	splineadvectionvx.cpp
	splitvlasovsolver.cpp
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
		DDC::DDC
		sll::splines
		vcx::speciesinfo
		vcx::geometry
)

add_library(vcx::vlasov ALIAS vlasov)
