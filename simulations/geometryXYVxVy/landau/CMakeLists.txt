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
            gyselalibxx::advection
            gyselalibxx::geometry_xyvxvy_${PRECISION}
            gyselalibxx::initialisation_xyvxvy_${PRECISION}
            gyselalibxx::interpolation
            gyselalibxx::io
            gyselalibxx::mpi_parallelisation
            gyselalibxx::paraconfpp
            gyselalibxx::poisson_xy_${PRECISION}
            gyselalibxx::quadrature
            gyselalibxx::time_integration_xyvxvy_${PRECISION}
            gyselalibxx::utils
            gyselalibxx::vlasov_xyvxvy_${PRECISION}
    )

    install(TARGETS ${target_name})
    string(TOUPPER ${INTERPOLATION} INTERPOLATION_MACRO)
    target_compile_definitions(${target_name} PRIVATE ${INTERPOLATION_MACRO})
    target_compile_definitions(${target_name} INTERFACE GYSELALIBXX_BUILD_REAL_PRECISION=${PRECISION})
endfunction()

build_landau4d_fft_target(spline double)
build_landau4d_fft_target(lagrange float)
build_landau4d_fft_target(lagrange double)
