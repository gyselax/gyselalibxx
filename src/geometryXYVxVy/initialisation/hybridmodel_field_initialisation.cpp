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

DFieldXY Hybridmodel_field_initialisation::operator()(DFieldXY const magnetic_field_z, DFieldXY const pressure) const
{
    IdxRangeXY const gridxy = get_idx_range<GridX, GridY>(magnetic_field_z);

    double const magnetic_kx = m_magnetic_init_perturb_mode * 2. * M_PI
                      / ddcHelper::total_interval_length(ddc::select<GridX>(gridxy));

    double const magnetic_ky = m_magnetic_init_perturb_mode * 2. * M_PI
                      / ddcHelper::total_interval_length(ddc::select<GridY>(gridxy));

    double const pressure_kx = m_pressure_init_perturb_mode * 2. * M_PI
                      / ddcHelper::total_interval_length(ddc::select<GridX>(gridxy));

    double const pressure_ky = m_pressure_init_perturb_mode * 2. * M_PI
                      / ddcHelper::total_interval_length(ddc::select<GridY>(gridxy));

    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            gridxy,
            KOKKOS_LAMBDA(IdxXY const ixy) {
                double const x = ddc::coordinate(ddc::select<GridX>(ixy));
                double const y = ddc::coordinate(ddc::select<GridY>(ixy));
                magnetic_field_z(ixy) = 1.0 + m_magnetic_init_perturb_amplitude * sin(magnetic_kx * x)
                                            + m_magnetic_init_perturb_amplitude * sin(magnetic_ky * y);
            });

    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            gridxy,
            KOKKOS_LAMBDA(IdxXY const ixy) {
                double const x = ddc::coordinate(ddc::select<GridX>(ixy));
                double const y = ddc::coordinate(ddc::select<GridY>(ixy));
                pressure(ixy) = 1.0 + m_pressure_init_perturb_amplitude * sin(pressure_kx * x)
                                    + m_pressure_init_perturb_amplitude * sin(pressure_ky * y);
            });

    return magnetic_field_z;
}
