add_executable(ddc_benchmarks deepcopy.cpp)
target_link_libraries(ddc_benchmarks
	PUBLIC
		benchmark::benchmark
		ddc::ddc
)
