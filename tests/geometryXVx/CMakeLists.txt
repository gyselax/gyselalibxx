# SPDX-License-Identifier: MIT

cmake_minimum_required(VERSION 3.15)

include(GoogleTest)

foreach(GEOMETRY_VARIANT IN LISTS GEOMETRY_XVx_VARIANTS_LIST)

add_executable(unit_tests_${GEOMETRY_VARIANT}
    collisions_inter.cpp
    collisions_intra_gridvx.cpp
    collisions_intra_maxwellian.cpp
    fluid_moments.cpp
    kineticsource.cpp
    krooksource.cpp
    ../main.cpp
    masks.cpp
    quadrature.cpp
    splitvlasovsolver.cpp
)
target_compile_features(unit_tests_${GEOMETRY_VARIANT} PUBLIC cxx_std_17)
target_link_libraries(unit_tests_${GEOMETRY_VARIANT}
    PUBLIC
        DDC::PDI_Wrapper
        GTest::gtest
        GTest::gmock
        paraconf::paraconf
        sll::splines
        vcx::advection
        vcx::boltzmann_${GEOMETRY_VARIANT}
        vcx::initialization_${GEOMETRY_VARIANT}
        vcx::poisson_${GEOMETRY_VARIANT}
        vcx::quadrature
        vcx::utils_${GEOMETRY_VARIANT}
)

gtest_discover_tests(unit_tests_${GEOMETRY_VARIANT}
    TEST_SUFFIX "_${GEOMETRY_VARIANT}")

endforeach()

target_sources(unit_tests_xperiod_vx
    PRIVATE
        femperiodicpoissonsolver.cpp
)

target_sources(unit_tests_xnonperiod_vx PRIVATE femnonperiodicpoissonsolver.cpp)

add_subdirectory(bump_on_tail)

add_subdirectory(landau)



foreach(GEOMETRY_VARIANT IN LISTS GEOMETRY_XVx_VARIANTS_LIST)
	add_executable(spline_quadrature_tests_${GEOMETRY_VARIANT}
	    ../main.cpp
	    quadrature_spline.cpp
	)
	target_compile_features(spline_quadrature_tests_${GEOMETRY_VARIANT} PUBLIC cxx_std_17)
	target_link_libraries(spline_quadrature_tests_${GEOMETRY_VARIANT}
	    PUBLIC
	        DDC::PDI_Wrapper
	        GTest::gtest
	        GTest::gmock
	        paraconf::paraconf
	        sll::splines
	        vcx::quadrature
	        vcx::geometry_${GEOMETRY_VARIANT}
	)
	
	gtest_discover_tests(spline_quadrature_tests_${GEOMETRY_VARIANT}
	    TEST_SUFFIX "_${GEOMETRY_VARIANT}")
endforeach()