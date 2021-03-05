cmake_minimum_required(VERSION 3.15)

add_executable(voicexx
	distributed2dfield.cpp
	simpleheat.cpp
)
target_include_directories(voicexx PRIVATE
	"${CMAKE_CURRENT_SOURCE_DIR}"
)
target_link_libraries(voicexx PRIVATE
		MPI::MPI_CXX
		std::mdspan
)
