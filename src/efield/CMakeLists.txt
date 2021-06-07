add_library(efield STATIC
	nullefieldsolver.cpp
)

target_compile_features(efield
	PUBLIC
		cxx_std_17
)

target_include_directories(efield
	PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries(efield
	PUBLIC
		ddc::ddc
		vcx::geometry
)

add_library(vcx::efield ALIAS efield)
