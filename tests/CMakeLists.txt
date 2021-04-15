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
