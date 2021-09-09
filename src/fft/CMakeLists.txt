add_library(fft INTERFACE
)

target_compile_features(fft
	INTERFACE
		cxx_std_17
)

target_include_directories(fft
	INTERFACE
		"${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries(fft
	INTERFACE
		ddc::ddc
                vcx::geometry
)


add_library(vcx::fft ALIAS fft)
