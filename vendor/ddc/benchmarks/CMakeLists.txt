# SPDX-License-Identifier: MIT

add_executable(ddc_benchmarks deepcopy.cpp)
target_link_libraries(ddc_benchmarks
	PUBLIC
		benchmark::benchmark
		DDC::DDC
)

add_executable(ddc_heat_equation_benchmark heat_equation.cpp)
target_link_libraries(ddc_heat_equation_benchmark
	PUBLIC
		benchmark::benchmark
		DDC::DDC
)
