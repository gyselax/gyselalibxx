cmake_minimum_required(VERSION 3.15)

add_subdirectory(array_tools)
add_subdirectory(boundary_conditions)
add_subdirectory(efield)
add_subdirectory(geometry)
add_subdirectory(math_tools)
add_subdirectory(matrix)
add_subdirectory(mesh)
add_subdirectory(splines)
add_subdirectory(time)
add_subdirectory(vlasov)

add_executable(voicexx
	voicexx.cpp
)

target_include_directories(voicexx
	PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries(voicexx
	PUBLIC
		vcx::array_tools
		sll::boundary_conditions
		vcx::efield
		vcx::mesh
		sll::splines
		vcx::time
		vcx::vlasov
)
