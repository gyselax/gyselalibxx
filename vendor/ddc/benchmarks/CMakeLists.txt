# SPDX-License-Identifier: MIT

add_executable(ddc_benchmark_deepcopy deepcopy.cpp)
target_link_libraries(ddc_benchmark_deepcopy
	PUBLIC
		benchmark::benchmark
		DDC::DDC
)

if("${BUILD_SPLINES_KERNEL}")
  add_executable(ddc_benchmark_splines splines.cpp)
  target_link_libraries(ddc_benchmark_splines
	  PUBLIC
		  benchmark::benchmark
		  DDC::DDC
  )
endif()
