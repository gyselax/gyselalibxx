# SPDX-License-Identifier: MIT
cmake_minimum_required(VERSION 3.15)
project(DDC VERSION 0.0.0 LANGUAGES CXX)

# List of options

option(BUILD_BENCHMARKS    "Build the benchmarks." OFF)
option(BUILD_DOCUMENTATION "Build the documentation" OFF)
option(DDC_BUILD_PDI_WRAPPER "Build the PDI wrapper" ON)

set(DDC_DEPENDENCY_POLICIES "AUTO" "EMBEDDED" "INSTALLED" "SUBPROJECT")

set(DDC_mdspan_DEPENDENCY_POLICY "AUTO" CACHE STRING "Policy to find the `mdspan` package. Options: ${DDC_DEPENDENCY_POLICIES}")
set_property(CACHE DDC_mdspan_DEPENDENCY_POLICY PROPERTY STRINGS "${DDC_DEPENDENCY_POLICIES}")

# Dependencies

if(DDC_mdspan_DEPENDENCY_POLICY STREQUAL "AUTO")
	if(NOT ${mdspan_FOUND})
		find_package(mdspan CONFIG)
	endif()
	if(NOT ${mdspan_FOUND})
		add_subdirectory(vendor/mdspan)
	endif()
elseif(DDC_mdspan_DEPENDENCY_POLICY STREQUAL "EMBEDDED")
	add_subdirectory(vendor/mdspan)
elseif(DDC_mdspan_DEPENDENCY_POLICY STREQUAL "INSTALLED")
	find_package(mdspan CONFIG REQUIRED)
endif()

# Our project

## The library itself
add_library(DDC INTERFACE)
target_compile_features(DDC INTERFACE cxx_std_17)
target_include_directories(DDC
	INTERFACE
		"$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
)
target_include_directories(DDC
	SYSTEM INTERFACE
		"$<INSTALL_INTERFACE:include>"
)
target_link_libraries(DDC
	INTERFACE std::mdspan)
add_library(DDC::DDC ALIAS DDC)
install(
	TARGETS DDC
	EXPORT DDCTargets)

## The PDI wrapper
if(DDC_BUILD_PDI_WRAPPER)
	if(NOT ${PDI_FOUND})
		find_package(PDI REQUIRED COMPONENTS C)
	endif()

	add_library(PDI INTERFACE)
	target_compile_features(PDI INTERFACE cxx_std_17)
	target_include_directories(PDI
		INTERFACE
			"$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
	)
	target_include_directories(PDI
		SYSTEM INTERFACE
			"$<INSTALL_INTERFACE:include>"
	)
	target_link_libraries(PDI
		INTERFACE
			DDC::DDC
			PDI::PDI_C)
	add_library(DDC::PDI ALIAS PDI)
	install(
		TARGETS PDI
		EXPORT DDCTargets)
endif()

include(CTest)

## if tests are enabled, build the tests in `tests/`
if(BUILD_TESTING)
	set(DDC_GTest_DEPENDENCY_POLICY "AUTO" CACHE STRING "Policy to find the `GTest` package. Options: ${DDC_DEPENDENCY_POLICIES}")
	set_property(CACHE DDC_GTest_DEPENDENCY_POLICY PROPERTY STRINGS ${DDC_DEPENDENCY_POLICIES})

	if(DDC_GTest_DEPENDENCY_POLICY STREQUAL "AUTO")
		if(NOT ${GTest_FOUND})
			find_package(GTest)
		endif()
		if(NOT ${GTest_FOUND})
			add_subdirectory(vendor/googletest)
		endif()
	elseif(DDC_GTest_DEPENDENCY_POLICY STREQUAL "EMBEDDED")
		add_subdirectory(vendor/googletest)
	elseif(DDC_GTest_DEPENDENCY_POLICY STREQUAL "INSTALLED")
		find_package(GTest REQUIRED)
	endif()

	add_subdirectory(tests/)
endif()

## if benchmarks are enabled, build the benchmarks in `benchmarks/`
if(BUILD_BENCHMARKS)
	set(DDC_benchmark_DEPENDENCY_POLICY "AUTO" CACHE STRING "Policy to find the `benchmark` package. Options: ${DDC_DEPENDENCY_POLICIES}")
	set_property(CACHE DDC_benchmark_DEPENDENCY_POLICY PROPERTY STRINGS ${DDC_DEPENDENCY_POLICIES})

	if(DDC_benchmark_DEPENDENCY_POLICY STREQUAL "AUTO")
		if(NOT ${benchmark_FOUND})
			find_package(benchmark)
		endif()
		if(NOT ${benchmark_FOUND})
			add_subdirectory(vendor/benchmark)
		endif()
	elseif(DDC_benchmark_DEPENDENCY_POLICY STREQUAL "EMBEDDED")
		add_subdirectory(vendor/benchmark)
	elseif(DDC_benchmark_DEPENDENCY_POLICY STREQUAL "INSTALLED")
		find_package(benchmark REQUIRED)
	endif()

	add_subdirectory(benchmarks/)
endif()

## if documentation is enabled, build the documentation docs
if(BUILD_DOCUMENTATION)
	add_subdirectory(docs/)
endif()

install(
	EXPORT DDCTargets
	NAMESPACE DDC::
	DESTINATION lib/cmake/DDC)

install(
	DIRECTORY include/
	TYPE INCLUDE)

include(CMakePackageConfigHelpers)
configure_package_config_file(cmake/DDCConfig.cmake.in
	${CMAKE_CURRENT_BINARY_DIR}/DDCConfig.cmake
	INSTALL_DESTINATION lib/cmake/DDC)

install(
	FILES ${CMAKE_CURRENT_BINARY_DIR}/DDCConfig.cmake
	DESTINATION lib/cmake/DDC)
