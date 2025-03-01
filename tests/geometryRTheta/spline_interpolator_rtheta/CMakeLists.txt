# SPDX-License-Identifier: MIT


set(test_name "spline_interpolator_rtheta_tests")


add_executable("${test_name}"
	spline_interpolator_rtheta.cpp
)
target_link_libraries("${test_name}"
    PUBLIC
        DDC::core
        gslx::geometry_RTheta
        gslx::interpolation_polar
        gslx::paraconfpp
        gslx::utils

)

 
find_package(Python3 REQUIRED COMPONENTS Interpreter)
 
add_test(NAME Test_spline_interpolator_rtheta_Convergence
       COMMAND "$<TARGET_FILE:Python3::Interpreter>" "${CMAKE_CURRENT_SOURCE_DIR}/test_convergence.py"
           "$<TARGET_FILE:${test_name}>")
           
set_property(TEST Test_spline_interpolator_rtheta_Convergence PROPERTY TIMEOUT 200)
set_property(TEST Test_spline_interpolator_rtheta_Convergence PROPERTY COST 100)
