add_library(efield STATIC
	efieldsolver.cpp
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
		vcx::mesh
)

add_library(vcx::efield ALIAS efield)
