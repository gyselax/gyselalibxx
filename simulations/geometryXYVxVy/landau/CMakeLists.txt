# SPDX-License-Identifier: MIT

set(INTERPOLATION_LIST spline lagrange)

foreach(INTERPOLATION IN LISTS INTERPOLATION_LIST)
    set(target_name landau4d_fft_${INTERPOLATION})
    add_executable(${target_name} landau4d_fft.cpp)
    target_link_libraries(${target_name}
        PUBLIC
            DDC::core
            DDC::pdi
            paraconf::paraconf
            PDI::PDI_C
            gslx::advection
            gslx::geometry_xyvxvy
            gslx::initialisation_xyvxvy
            gslx::interpolation
            gslx::io
            gslx::mpi_parallelisation
            gslx::paraconfpp
            gslx::poisson_xy
            gslx::quadrature
            gslx::time_integration_xyvxvy
            gslx::utils
            gslx::vlasov_xyvxvy
    )

    install(TARGETS ${target_name})
    string(TOUPPER ${INTERPOLATION} INTERPOLATION_MACRO)
    target_compile_definitions(${target_name} PRIVATE ${INTERPOLATION_MACRO})
endforeach()
