# SPDX-License-Identifier: MIT

add_executable(unit_tests_quadrature
    quadrature_1d.cpp
    quadrature_2d.cpp
	neumann_quadrature_spline.cpp
    batched_quadrature.cpp
    ../main.cpp
)
target_link_libraries(unit_tests_quadrature
    PUBLIC
        DDC::DDC
        GTest::gtest
        GTest::gmock

        sll::SLL

        gslx::quadrature
        gslx::utils
)
gtest_discover_tests(unit_tests_quadrature 
		DISCOVERY_MODE PRE_TEST)
target_compile_definitions(unit_tests_quadrature PUBLIC)

foreach(GEOMETRY_VARIANT IN LISTS GEOMETRY_XVx_VARIANTS_LIST)
	add_executable(spline_quadrature_tests_${GEOMETRY_VARIANT}
	    ../main.cpp
	    quadrature_spline.cpp
	)
	target_link_libraries(spline_quadrature_tests_${GEOMETRY_VARIANT}
	    PUBLIC
            DDC::DDC
	        GTest::gtest
	        GTest::gmock
	        paraconf::paraconf
	        sll::SLL
	        gslx::quadrature
            gslx::utils
	)

	gtest_discover_tests(spline_quadrature_tests_${GEOMETRY_VARIANT}
	    TEST_SUFFIX "_${GEOMETRY_VARIANT}"
        DISCOVERY_MODE PRE_TEST)
endforeach()
target_compile_definitions(spline_quadrature_tests_xperiod_vx INTERFACE PERIODIC_RDIMX)
