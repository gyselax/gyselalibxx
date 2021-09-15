cmake_minimum_required(VERSION 3.15)
project(ddc CXX)

# List of options
option(BUILD_BENCHMARKS "Build the benchmarks." OFF)

# Our dependencies

## Use CTest for running tests
include(CTest)

# Our project

## The library itself
add_library(ddc INTERFACE)
target_compile_features(ddc INTERFACE cxx_std_17)
target_include_directories(ddc
	INTERFACE
		"$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
)
target_link_libraries(ddc
	INTERFACE std::mdspan)
add_library(ddc::ddc ALIAS ddc)

## The PDI wrapper
add_library(ddc_pdi INTERFACE)
target_compile_features(ddc_pdi INTERFACE cxx_std_17)
target_include_directories(ddc_pdi
	INTERFACE
		"$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
)
target_link_libraries(ddc_pdi
	INTERFACE
		std::mdspan
		PDI::PDI_C)
add_library(ddc::PDI ALIAS ddc_pdi)

## if tests are enabled, build the tests in `tests/`
if("${BUILD_TESTING}")
    add_subdirectory(tests/)
endif()

## if benchmarks are enabled, build the benchmarks in `benchmarks/`
if("${BUILD_BENCHMARKS}")
    add_subdirectory(benchmarks/)
endif()
