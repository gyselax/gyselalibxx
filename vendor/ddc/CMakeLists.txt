# Copyright (C) The DDC development team, see COPYRIGHT.md file
#
# SPDX-License-Identifier: MIT

cmake_minimum_required(VERSION 3.15)
project(DDC VERSION 0.0.0 LANGUAGES CXX)

# List of options

option(DDC_BUILD_BENCHMARKS       "Build DDC benchmarks." OFF)
option(DDC_BUILD_DOCUMENTATION    "Build DDC documentation/website" OFF)
option(DDC_BUILD_DOUBLE_PRECISION "Build DDC with double precision support, float is used otherwise" ON)
option(DDC_BUILD_EXAMPLES         "Build DDC examples" ON)
option(DDC_BUILD_KERNELS_FFT      "Build DDC kernels for FFT" ON)
option(DDC_BUILD_KERNELS_SPLINES  "Build DDC kernels for splines" ON)
option(DDC_BUILD_PDI_WRAPPER      "Build DDC PDI wrapper" ON)
option(DDC_BUILD_TESTS            "Build DDC tests if BUILD_TESTING is enabled" ON)

# Dependencies

set(DDC_DEPENDENCY_POLICIES "AUTO" "EMBEDDED" "INSTALLED" "SUBPROJECT")


## CMake modules

include(CMakePackageConfigHelpers)
include(CTest)


## mdspan

set(DDC_mdspan_DEPENDENCY_POLICY "AUTO" CACHE STRING "Policy to find the `mdspan` package. Options: ${DDC_DEPENDENCY_POLICIES}")
set_property(CACHE DDC_mdspan_DEPENDENCY_POLICY PROPERTY STRINGS "${DDC_DEPENDENCY_POLICIES}")
if("${DDC_mdspan_DEPENDENCY_POLICY}" STREQUAL "AUTO")
	if(NOT TARGET std::mdspan)
		find_package(mdspan 0.6.0 EXACT CONFIG QUIET)
		if(NOT mdspan_FOUND)
			add_subdirectory(vendor/mdspan)
		endif()
	endif()
elseif("${DDC_mdspan_DEPENDENCY_POLICY}" STREQUAL "EMBEDDED")
	add_subdirectory(vendor/mdspan)
elseif("${DDC_mdspan_DEPENDENCY_POLICY}" STREQUAL "INSTALLED")
	find_package(mdspan 0.6.0 EXACT CONFIG REQUIRED)
endif()

## kokkos

set(DDC_Kokkos_DEPENDENCY_POLICY "AUTO" CACHE STRING "Policy to find the `Kokkos` package. Options: ${DDC_DEPENDENCY_POLICIES}")
set_property(CACHE DDC_Kokkos_DEPENDENCY_POLICY PROPERTY STRINGS "${DDC_DEPENDENCY_POLICIES}")
if("${DDC_Kokkos_DEPENDENCY_POLICY}" STREQUAL "AUTO")
	if(NOT TARGET Kokkos::kokkos)
		find_package(Kokkos CONFIG QUIET)
		if(NOT Kokkos_FOUND)
			if("${Kokkos_ENABLE_CUDA}")
				option(Kokkos_ENABLE_CUDA_CONSTEXPR "Whether to activate experimental relaxed constexpr functions" ON)
				if(NOT "${Kokkos_ENABLE_CUDA_CONSTEXPR}")
					message(SEND_ERROR "Kokkos_ENABLE_CUDA_CONSTEXPR must be enabled, currently it is not")
				endif()
				option(Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE "Whether to enable relocatable device code (RDC) for CUDA" ON)
				if(NOT "${Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE}")
					message(SEND_ERROR "Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE must be enabled, currently it is not")
				endif()
			endif()
			if("${Kokkos_ENABLE_HIP}")
				option(Kokkos_ENABLE_HIP_RELOCATABLE_DEVICE_CODE "Whether to enable relocatable device code (RDC) for HIP" ON)
				if(NOT "${Kokkos_ENABLE_HIP_RELOCATABLE_DEVICE_CODE}")
					message(SEND_ERROR "Kokkos_ENABLE_HIP_RELOCATABLE_DEVICE_CODE must be enabled, currently it is not")
				endif()
			endif()
			add_subdirectory(vendor/kokkos)
		endif()
	endif()
elseif("${DDC_Kokkos_DEPENDENCY_POLICY}" STREQUAL "EMBEDDED")
	if("${Kokkos_ENABLE_CUDA}")
		option(Kokkos_ENABLE_CUDA_CONSTEXPR "Whether to activate experimental relaxed constexpr functions" ON)
		if(NOT "${Kokkos_ENABLE_CUDA_CONSTEXPR}")
			message(SEND_ERROR "Kokkos_ENABLE_CUDA_CONSTEXPR must be enabled, currently it is not")
		endif()
		option(Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE "Whether to enable relocatable device code (RDC) for CUDA" ON)
		if(NOT "${Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE}")
			message(SEND_ERROR "Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE must be enabled, currently it is not")
		endif()
	endif()
	if("${Kokkos_ENABLE_HIP}")
		option(Kokkos_ENABLE_HIP_RELOCATABLE_DEVICE_CODE "Whether to enable relocatable device code (RDC) for HIP" ON)
		if(NOT "${Kokkos_ENABLE_HIP_RELOCATABLE_DEVICE_CODE}")
			message(SEND_ERROR "Kokkos_ENABLE_HIP_RELOCATABLE_DEVICE_CODE must be enabled, currently it is not")
		endif()
	endif()
	add_subdirectory(vendor/kokkos)
elseif("${DDC_Kokkos_DEPENDENCY_POLICY}" STREQUAL "INSTALLED")
	find_package(Kokkos CONFIG REQUIRED)
endif()

# FFTW
list( APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" ) # Maybe not specific to FFTW
if("${DDC_BUILD_KERNELS_FFT}" AND NOT FFTW_FOUND)
	find_package( FFTW MODULE REQUIRED )
endif()

## CUDA + CUDAToolkit
if("${DDC_BUILD_KERNELS_FFT}" AND "${Kokkos_ENABLE_CUDA}")
	find_package( CUDAToolkit MODULE REQUIRED )
	if( NOT(CUDAToolkit_FOUND) )
	  message(FATAL_ERROR "CUDAToolkit not found." )
	endif()
endif()


## PDI

if("${DDC_BUILD_PDI_WRAPPER}")
	if(NOT TARGET PDI::PDI_C)
		find_package(PDI REQUIRED COMPONENTS C)
	endif()
endif()


## GoogleTest

if("${BUILD_TESTING}" AND "${DDC_BUILD_TESTS}")
	set(DDC_GTest_DEPENDENCY_POLICY "AUTO" CACHE STRING "Policy to find the `GTest` package. Options: ${DDC_DEPENDENCY_POLICIES}")
	set_property(CACHE DDC_GTest_DEPENDENCY_POLICY PROPERTY STRINGS ${DDC_DEPENDENCY_POLICIES})
	set(GOOGLETEST_PATH "${CMAKE_CURRENT_SOURCE_DIR}/vendor/googletest" CACHE PATH "Path to the googletest library source")

	if("${DDC_GTest_DEPENDENCY_POLICY}" STREQUAL "AUTO")
		if(NOT TARGET GTest::GTest AND NOT TARGET GTest::gtest)
			find_package(GTest)
			if(NOT GTest_FOUND)
				add_subdirectory("${GOOGLETEST_PATH}")
			endif()
		endif()
		if(NOT TARGET GTest::gtest_main)
      if(TARGET GTest::GTest)
			  add_library(GTest::gtest_main ALIAS GTest::GTest)
      else()
			  add_library(GTest::gtest_main ALIAS GTest::gtest)
      endif()
		endif()
	elseif("${DDC_GTest_DEPENDENCY_POLICY}" STREQUAL "EMBEDDED")
		add_subdirectory("${GOOGLETEST_PATH}")
	elseif("${DDC_GTest_DEPENDENCY_POLICY}" STREQUAL "INSTALLED")
		find_package(GTest REQUIRED)
		if(NOT TARGET GTest::gtest_main)
			add_library(GTest::gtest_main ALIAS GTest::GTest)
		endif()
	endif()

	include(GoogleTest)
endif()


## Google Benchmark

if("${DDC_BUILD_BENCHMARKS}")
	set(DDC_benchmark_DEPENDENCY_POLICY "AUTO" CACHE STRING "Policy to find the `benchmark` package. Options: ${DDC_DEPENDENCY_POLICIES}")
	set_property(CACHE DDC_benchmark_DEPENDENCY_POLICY PROPERTY STRINGS ${DDC_DEPENDENCY_POLICIES})
	option(BENCHMARK_ENABLE_TESTING "Enable testing of the benchmark library." OFF)
	option(BENCHMARK_ENABLE_INSTALL "Enable installation of benchmark. (Projects embedding benchmark may want to turn this OFF.)" OFF)

	if("${DDC_benchmark_DEPENDENCY_POLICY}" STREQUAL "AUTO")
		if(NOT TARGET benchmark::benchmark)
			find_package(benchmark QUIET)
			if(NOT benchmark_FOUND)
				add_subdirectory(vendor/benchmark)
			endif()
		endif()
	elseif("${DDC_benchmark_DEPENDENCY_POLICY}" STREQUAL "EMBEDDED")
		add_subdirectory(vendor/benchmark)
	elseif("${DDC_benchmark_DEPENDENCY_POLICY}" STREQUAL "INSTALLED")
		find_package(benchmark REQUIRED)
	endif()
endif()


## Doxygen

if("${DDC_BUILD_DOCUMENTATION}")
	find_package(Doxygen 1.8.13 REQUIRED OPTIONAL_COMPONENTS dot)
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
	INTERFACE std::mdspan Kokkos::kokkos
)
target_compile_definitions(DDC
	INTERFACE
		MDSPAN_USE_PAREN_OPERATOR=1
)
if("${DDC_BUILD_DOUBLE_PRECISION}")
	target_compile_definitions(DDC INTERFACE DDC_BUILD_DOUBLE_PRECISION)
endif()
add_library(DDC::DDC ALIAS DDC)
install(
	TARGETS DDC
	EXPORT DDCTargets
)

# Link library to DDC

if( FFTW_FOUND )
	target_link_libraries( DDC INTERFACE FFTW::Float )
	target_link_libraries( DDC INTERFACE FFTW::Double )
endif()

if( FFTW_FOUND AND "${Kokkos_ENABLE_SERIAL}")
	target_compile_definitions(DDC INTERFACE fftw_serial_AVAIL)
endif()

if( FFTW_FOUND AND "${Kokkos_ENABLE_OPENMP}")
	target_link_libraries( DDC INTERFACE FFTW::FloatOpenMP )
	target_link_libraries( DDC INTERFACE FFTW::DoubleOpenMP )
	target_compile_definitions(DDC INTERFACE fftw_omp_AVAIL)
endif()

if( CUDAToolkit_FOUND )
	target_link_libraries( DDC INTERFACE CUDA::cufft )
	target_compile_definitions(DDC INTERFACE cufft_AVAIL)
endif()

if("${Kokkos_ENABLE_HIP}")
	find_package( hipfft REQUIRED ) # Not compatible with nvidia
	target_link_libraries( DDC INTERFACE hip::hipfft )
	target_compile_definitions(DDC INTERFACE hipfft_AVAIL)
endif()

if("${DDC_BUILD_KERNELS_SPLINES}")
  # Ginkgo
  find_package(Ginkgo 1.7.0 EXACT REQUIRED)
  target_link_libraries(DDC INTERFACE Ginkgo::ginkgo)
  target_compile_definitions(DDC INTERFACE ginkgo_AVAIL)
endif()

## The PDI wrapper

if("${DDC_BUILD_PDI_WRAPPER}")
	add_library(PDI_Wrapper INTERFACE)
	target_compile_features(PDI_Wrapper INTERFACE cxx_std_17)
	target_include_directories(PDI_Wrapper
		INTERFACE
			"$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
	)
	target_include_directories(PDI_Wrapper
		SYSTEM INTERFACE
			"$<INSTALL_INTERFACE:include>"
	)
	target_link_libraries(PDI_Wrapper
		INTERFACE
			DDC::DDC
			PDI::PDI_C)
	target_compile_definitions(PDI_Wrapper INTERFACE "DDC_BUILD_PDI_WRAPPER")
	add_library(DDC::PDI_Wrapper ALIAS PDI_Wrapper)
	install(
		TARGETS PDI_Wrapper
		EXPORT DDCTargets)
endif()


## if examples are enabled, build them

if("${DDC_BUILD_EXAMPLES}")
	add_subdirectory(examples/)
endif()


## if tests are enabled, build them

if("${BUILD_TESTING}" AND "${DDC_BUILD_TESTS}")
	add_subdirectory(tests/)
endif()

## if benchmarks are enabled, build them

if("${DDC_BUILD_BENCHMARKS}")
	add_subdirectory(benchmarks/)
endif()


## if documentation is enabled, build it

if("${DDC_BUILD_DOCUMENTATION}")
	add_subdirectory(docs/)
endif()


## installation

install(
	EXPORT DDCTargets
	NAMESPACE DDC::
	DESTINATION lib/cmake/DDC)

install(
	DIRECTORY include/
	TYPE INCLUDE)

configure_package_config_file(cmake/DDCConfig.cmake.in
	${CMAKE_CURRENT_BINARY_DIR}/DDCConfig.cmake
	INSTALL_DESTINATION lib/cmake/DDC)

install(
	FILES ${CMAKE_CURRENT_BINARY_DIR}/DDCConfig.cmake
	DESTINATION lib/cmake/DDC)
