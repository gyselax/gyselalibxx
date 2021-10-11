cmake_minimum_required(VERSION 3.15)

add_subdirectory(efield)
add_subdirectory(fft)
add_subdirectory(geometry)
add_subdirectory(time)
add_subdirectory(fdistribu)
add_subdirectory(vlasov)

add_executable(voicexx
	voicexx.cpp
)
target_compile_features(voicexx PUBLIC cxx_std_17)
target_link_libraries(voicexx
	PUBLIC
		DDC::DDC
		DDC::PDI
		vcx::fft
		vcx::efield
		vcx::time
		vcx::fdistribu
		vcx::vlasov
		FFTW::Double
		PDI::pdi
)
