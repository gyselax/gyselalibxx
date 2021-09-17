add_library(deprecated_fft STATIC
	fft_1d.cpp
)

target_compile_features(deprecated_fft
	PUBLIC
		cxx_std_17
)

target_include_directories(deprecated_fft
	PUBLIC
		"${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries(deprecated_fft
	PUBLIC
		DDC::DDC
		vcx::geometry
)


add_library(vcx::deprecated::fft ALIAS deprecated_fft)
