add_library(fft STATIC
        fft_1d.cpp
)

target_compile_features(fft
	PUBLIC
		cxx_std_17
)

target_include_directories(fft
	PUBLIC
		"${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries(fft
	PUBLIC
		ddc::ddc
                vcx::geometry
)


add_library(vcx::fft ALIAS fft)
