// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include "hybridmodel_field_initialisation.hpp"

Hybridmodel_field_initialisation::Hybridmodel_field_initialisation(
        int magnetic_init_perturb_mode,
        double magnetic_init_perturb_amplitude,
        int pressure_init_perturb_mode,
        double pressure_init_perturb_amplitude)
    : m_magnetic_init_perturb_mode(std::move(magnetic_init_perturb_mode))
    , m_magnetic_init_perturb_amplitude(std::move(magnetic_init_perturb_amplitude))
    , m_pressure_init_perturb_mode(std::move(pressure_init_perturb_mode))
    , m_pressure_init_perturb_amplitude(std::move(pressure_init_perturb_amplitude))
{
}

DFieldX Hybridmodel_field_initialisation::operator()(DFieldX const magnetic_field_x, DFieldX const magnetic_field_y, DFieldX const magnetic_field_z, DFieldX const pressure) const
{
    IdxRangeX const gridx = get_idx_range<GridX>(magnetic_field_z);

    double const magnetic_kx = m_magnetic_init_perturb_mode * 2. * M_PI
                      / ddcHelper::total_interval_length(ddc::select<GridX>(gridx));


    double const pressure_kx = m_pressure_init_perturb_mode * 2. * M_PI
                      / ddcHelper::total_interval_length(ddc::select<GridX>(gridx));

    

    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            gridx,
            KOKKOS_LAMBDA(IdxX const ix) {
                double const x = ddc::coordinate(ddc::select<GridX>(ix));
                magnetic_field_x(ix) = 0.6;
                magnetic_field_y(ix) = 0.7 + m_magnetic_init_perturb_amplitude * cos(magnetic_kx * x); 
                magnetic_field_z(ix) = 1.0 + m_magnetic_init_perturb_amplitude * sin(magnetic_kx * x);  
                //magnetic_field_x(ix) = 0.6;
                //magnetic_field_y(ix) = 0.7;// + m_magnetic_init_perturb_amplitude * cos(magnetic_kx * x); 
                //magnetic_field_z(ix) = 1.0;// + m_magnetic_init_perturb_amplitude * sin(magnetic_kx * x);      
            });

    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            gridx,
            KOKKOS_LAMBDA(IdxX const ix) {
                double const x = ddc::coordinate(ddc::select<GridX>(ix));
                pressure(ix) = 1.0 + m_pressure_init_perturb_amplitude * sin(pressure_kx * x);
            });

    return magnetic_field_z;
}
