add_library(vlasov STATIC
	bsl_advection_x.cpp
	bsl_advection_vx.cpp
	nulladvectionvx.cpp
	splitvlasovsolver.cpp
	spline_interpolator_x.cpp
	spline_interpolator_vx.cpp
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
