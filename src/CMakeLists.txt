cmake_minimum_required(VERSION 3.15)

add_subdirectory(deprecated)
add_subdirectory(efield)
add_subdirectory(fft)
add_subdirectory(geometry)
add_subdirectory(time)
add_subdirectory(vlasov)

add_executable(voicexx
	voicexx.cpp
	expose_to_pdi.cpp
)
target_compile_features(voicexx PUBLIC cxx_std_17)
target_link_libraries(voicexx
	PUBLIC
		ddc::ddc
		vcx::fft
		vcx::efield
		vcx::time
		vcx::vlasov
		FFTW::Double
		PDI::pdi
)
