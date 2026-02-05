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

    double B_perturbation = 0.001;

    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            gridxy,
            KOKKOS_LAMBDA(IdxXY const ixy) {
                double const x = ddc::coordinate(ddc::select<GridX>(ixy));
                double const y = ddc::coordinate(ddc::select<GridY>(ixy));
                magnetic_field_z(ixy) = 1.0 + 0.0 * m_magnetic_init_perturb_amplitude * sin(magnetic_kx * x)
                                            + 0.0 * m_magnetic_init_perturb_amplitude * sin(magnetic_ky * y);
                
                magnetic_field_z(ixy) +=     B_perturbation * sin(magnetic_kx * x);  
                magnetic_field_z(ixy) +=     B_perturbation * sin(2.0 * magnetic_kx * x);  
                magnetic_field_z(ixy) +=     B_perturbation * sin(3.0 * magnetic_kx * x);  
                magnetic_field_z(ixy) +=     B_perturbation * sin(4.0 * magnetic_kx * x);  
                magnetic_field_z(ixy) +=     B_perturbation * sin(5.0 * magnetic_kx * x);  
                magnetic_field_z(ixy) +=     B_perturbation * sin(6.0 * magnetic_kx * x);  
                magnetic_field_z(ixy) +=     B_perturbation * sin(7.0 * magnetic_kx * x);  
                magnetic_field_z(ixy) +=     B_perturbation * sin(8.0 * magnetic_kx * x);  
                magnetic_field_z(ixy) +=     B_perturbation * sin(9.0 * magnetic_kx * x);  
                magnetic_field_z(ixy) +=     B_perturbation * sin(10.0 * magnetic_kx * x);  
                magnetic_field_z(ixy) +=     B_perturbation * sin(11.0 * magnetic_kx * x);  
                magnetic_field_z(ixy) +=     B_perturbation * sin(12.0 * magnetic_kx * x);  
                magnetic_field_z(ixy) +=     B_perturbation * sin(13.0 * magnetic_kx * x);  
                magnetic_field_z(ixy) +=     B_perturbation * sin(14.0 * magnetic_kx * x);  
                magnetic_field_z(ixy) +=     B_perturbation * sin(15.0 * magnetic_kx * x);  
                magnetic_field_z(ixy) +=     B_perturbation * sin(16.0 * magnetic_kx * x);  
                magnetic_field_z(ixy) +=     B_perturbation * sin(17.0 * magnetic_kx * x);  
                magnetic_field_z(ixy) +=     B_perturbation * sin(18.0 * magnetic_kx * x);  
                magnetic_field_z(ixy) +=     B_perturbation * sin(19.0 * magnetic_kx * x);  
                magnetic_field_z(ixy) +=     B_perturbation * sin(20.0 * magnetic_kx * x);  
                magnetic_field_z(ixy) +=     B_perturbation * sin(21.0 * magnetic_kx * x);  
                magnetic_field_z(ixy) +=     B_perturbation * sin(22.0 * magnetic_kx * x);  
                magnetic_field_z(ixy) +=     B_perturbation * sin(23.0 * magnetic_kx * x);  
                magnetic_field_z(ixy) +=     B_perturbation * sin(24.0 * magnetic_kx * x);  
                magnetic_field_z(ixy) +=     B_perturbation * sin(25.0 * magnetic_kx * x);  
                magnetic_field_z(ixy) +=     B_perturbation * sin(26.0 * magnetic_kx * x);  
                magnetic_field_z(ixy) +=     B_perturbation * sin(27.0 * magnetic_kx * x);  
                magnetic_field_z(ixy) +=     B_perturbation * sin(28.0 * magnetic_kx * x);  
                magnetic_field_z(ixy) +=     B_perturbation * sin(29.0 * magnetic_kx * x);  
                magnetic_field_z(ixy) +=     B_perturbation * sin(30.0 * magnetic_kx * x);  
                magnetic_field_z(ixy) +=     B_perturbation * sin(31.0 * magnetic_kx * x);  
                magnetic_field_z(ixy) +=     B_perturbation * sin(32.0 * magnetic_kx * x);  
            });

    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            gridxy,
            KOKKOS_LAMBDA(IdxXY const ixy) {
                double const x = ddc::coordinate(ddc::select<GridX>(ixy));
                double const y = ddc::coordinate(ddc::select<GridY>(ixy));
                pressure(ixy) = 0.09 + 0.0 * m_pressure_init_perturb_amplitude * sin(pressure_kx * x)
                                    + 0.0 * m_pressure_init_perturb_amplitude * sin(pressure_ky * y);
            });

    return magnetic_field_z;
}
