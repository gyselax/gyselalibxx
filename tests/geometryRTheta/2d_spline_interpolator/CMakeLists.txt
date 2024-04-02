# SPDX-License-Identifier: MIT

cmake_minimum_required(VERSION 3.15)

set(test_name "2d_spline_interpolator_tests")


add_executable(2d_spline_interpolator_tests
	2d_spline_interpolator.cpp
)
target_compile_features(2d_spline_interpolator_tests PUBLIC cxx_std_17)
target_link_libraries(2d_spline_interpolator_tests
    PUBLIC
        DDC::DDC
        sll::SLL
        gslx::geometry_RTheta
        gslx::interpolation_2D_rp
        gslx::paraconfpp
)

 
find_package(Python3 REQUIRED COMPONENTS Interpreter)
 
add_test(NAME Test_2d_spline_interpolator_Convergence
       COMMAND "$<TARGET_FILE:Python3::Interpreter>" "${CMAKE_CURRENT_SOURCE_DIR}/test_convergence.py"
           "$<TARGET_FILE:${test_name}>")
           
set_property(TEST Test_2d_spline_interpolator_Convergence PROPERTY TIMEOUT 200)
set_property(TEST Test_2d_spline_interpolator_Convergence PROPERTY COST 100)
