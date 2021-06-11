cmake_minimum_required(VERSION 3.15)

add_library(ddc STATIC
	src/deprecated/bsplines.cpp
	src/deprecated/bsplines_uniform.cpp
	src/deprecated/bsplines_non_uniform.cpp
	src/deprecated/spline_1d.cpp
	src/deprecated/spline_2d.cpp
	src/null_boundary_value.cpp
)

target_compile_features(ddc
	PUBLIC
		cxx_std_17
)

target_include_directories(ddc
	PUBLIC
		"$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
)

target_link_libraries(ddc
	PUBLIC
		mdspan_wrapper::mdsw
                sll::matrix
)

add_library(ddc::ddc ALIAS ddc)

add_subdirectory(tests/)

if(${BUILD_BENCHMARKS})
    add_subdirectory(benchmarks)
endif()