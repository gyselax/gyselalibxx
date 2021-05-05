cmake_minimum_required(VERSION 3.15)

add_executable(tst
	tst.cpp
)

target_compile_features(tst
	PUBLIC
		cxx_std_17
)

target_link_libraries(tst
	PUBLIC
		vcx::array_tools
)

add_test(NAME tst COMMAND tst)


add_executable(test_splines
	test_splines.cpp
)

target_compile_features(test_splines
	PUBLIC
		cxx_std_17
)

target_include_directories(test_splines
	PUBLIC
		"${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries(test_splines
	PUBLIC
                sll::boundary_conditions
)

add_test(NAME test_splines COMMAND test_splines)
