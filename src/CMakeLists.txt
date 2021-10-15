cmake_minimum_required(VERSION 3.15)

add_subdirectory(poisson)
add_subdirectory(fft)
add_subdirectory(geometry)
add_subdirectory(time)
add_subdirectory(speciesinfo)
add_subdirectory(vlasov)
add_subdirectory(initialization)

add_executable(voicexx
	voicexx.cpp
)
target_compile_features(voicexx PUBLIC cxx_std_17)
target_link_libraries(voicexx
	PUBLIC
		DDC::DDC
		DDC::PDI
		vcx::initialization
		vcx::fft
		vcx::poisson
		vcx::time
		vcx::speciesinfo
		vcx::vlasov
		FFTW::Double
		PDI::pdi
		sll::splines
)
