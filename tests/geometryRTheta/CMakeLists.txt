# SPDX-License-Identifier: MIT

# Allow the possibility of deactivating slow polar poisson tests
option(POISSON_2D_BUILD_TESTING "Build polar poisson tests" ON)


add_subdirectory(advection_field_rtheta)
add_subdirectory(advection_rtheta)
if("${POISSON_2D_BUILD_TESTING}")
    add_subdirectory(polar_poisson)
endif()
add_subdirectory(quadrature)
add_subdirectory(spline_interpolator_rtheta)
add_subdirectory(time_solver)
