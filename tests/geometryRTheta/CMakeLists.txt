# SPDX-License-Identifier: MIT

# Allow the possibility of deactivating slow polar poisson tests
option(POISSON_2D_BUILD_TESTING "Build polar poisson tests" ON)

if("${POISSON_2D_BUILD_TESTING}")
    add_subdirectory(polar_poisson)
endif()
add_subdirectory(2d_spline_interpolator)
add_subdirectory(quadrature)
add_subdirectory(advection_2d_rp)
add_subdirectory(advection_field_rp)
