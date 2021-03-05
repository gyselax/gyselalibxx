cmake_minimum_required(VERSION 3.15)

add_executable(voicexx
	cartesiandistribution2d.cpp
	commandlineconfig.cpp
	distributed2dfield.cpp
	finitediffheatsolver.cpp
	fixedconditions.cpp
	screenprinter.cpp
	simpleheat.cpp
	simulation.cpp
)
target_include_directories(voicexx PRIVATE
	"${CMAKE_CURRENT_SOURCE_DIR}"
)
target_link_libraries(voicexx PRIVATE
		MPI::MPI_CXX
		std::mdspan
)
