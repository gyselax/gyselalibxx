# SPDX-License-Identifier: MIT

function (build_landau4d_fft_target INTERPOLATION PRECISION)
    set(target_name landau4d_fft_${INTERPOLATION}_${PRECISION})
    add_executable(${target_name} landau4d_fft.cpp)
    target_link_libraries(${target_name}
        PUBLIC
            DDC::core
            DDC::pdi
            paraconf::paraconf
            PDI::PDI_C
            gslx::advection
            gslx::geometry_xyvxvy_${PRECISION}
            gslx::initialisation_xyvxvy_${PRECISION}
            gslx::interpolation
            gslx::io
            gslx::mpi_parallelisation
            gslx::paraconfpp
            gslx::poisson_xy_${PRECISION}
            gslx::quadrature
            gslx::time_integration_xyvxvy_${PRECISION}
            gslx::utils
            gslx::vlasov_xyvxvy_${PRECISION}
    )

    install(TARGETS ${target_name})
    string(TOUPPER ${INTERPOLATION} INTERPOLATION_MACRO)
    target_compile_definitions(${target_name} PRIVATE ${INTERPOLATION_MACRO})
    target_compile_definitions(${target_name} INTERFACE GYSELALIBXX_BUILD_REAL_PRECISION=${PRECISION})
endfunction()

build_landau4d_fft_target(spline double)
build_landau4d_fft_target(lagrange float)
build_landau4d_fft_target(lagrange double)
