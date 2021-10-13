add_library(poisson STATIC
	nullpoissonsolver.cpp
	fftpoissonsolver.cpp
)

target_compile_features(poisson
	PUBLIC
		cxx_std_17
)

target_include_directories(poisson
	PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries(poisson
	PUBLIC
		DDC::DDC
		sll::splines
		vcx::speciesinfo
		vcx::fft
		vcx::geometry
)

add_library(vcx::poisson ALIAS poisson)
