add_library(boundary_conditions STATIC
	null_boundary_value.cpp
)

target_include_directories(boundary_conditions
	PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

add_library(sll::boundary_conditions ALIAS boundary_conditions)
