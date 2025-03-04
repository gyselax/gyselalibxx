# SPDX-License-Identifier: MIT

# Allow the possibility of deactivating slow polar poisson tests
option(POISSON_2D_BUILD_TESTING "Build polar poisson tests" ON)

if("${POISSON_2D_BUILD_TESTING}")
    add_subdirectory(polar_poisson)
endif()
add_subdirectory(spline_interpolator_rtheta)
add_subdirectory(quadrature)
add_subdirectory(advection_rtheta)
add_subdirectory(advection_field_rtheta)
