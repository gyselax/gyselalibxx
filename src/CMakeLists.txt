cmake_minimum_required(VERSION 3.15)

add_subdirectory(efield)
add_subdirectory(fft)
add_subdirectory(geometry)
add_subdirectory(splines)
add_subdirectory(time)
add_subdirectory(vlasov)

add_executable(voicexx
	voicexx.cpp
)
target_compile_features(voicexx PUBLIC cxx_std_17)
target_link_libraries(voicexx
	PUBLIC
		ddc::ddc
		vcx::fft
		sll::splines
		vcx::efield
		vcx::time
		vcx::vlasov
		FFTW::Double
)
