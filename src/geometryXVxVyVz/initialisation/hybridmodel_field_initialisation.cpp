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

    
    //double B_perturbation = 0.001;

    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            gridx,
            KOKKOS_LAMBDA(IdxX const ix) {
                double const x = ddc::coordinate(ddc::select<GridX>(ix));
                
                magnetic_field_x(ix) = 1.0 + 0.0 * m_magnetic_init_perturb_amplitude *  Kokkos::sin(magnetic_kx * x) *  Kokkos::cos(magnetic_kx * x); ;
                magnetic_field_y(ix) = 0. - Kokkos::sin(4.0 * magnetic_kx * x); 
                magnetic_field_z(ix) = 0. + Kokkos::cos(4.0 * magnetic_kx * x);
                /*
                magnetic_field_z(ix) +=     B_perturbation * sin(2.0 * magnetic_kx * x);  
                magnetic_field_z(ix) +=     B_perturbation * sin(3.0 * magnetic_kx * x);  
                magnetic_field_z(ix) +=     B_perturbation * sin(4.0 * magnetic_kx * x);  
                magnetic_field_z(ix) +=     B_perturbation * sin(5.0 * magnetic_kx * x);  
                magnetic_field_z(ix) +=     B_perturbation * sin(6.0 * magnetic_kx * x);  
                magnetic_field_z(ix) +=     B_perturbation * sin(7.0 * magnetic_kx * x);  
                magnetic_field_z(ix) +=     B_perturbation * sin(8.0 * magnetic_kx * x);  
                magnetic_field_z(ix) +=     B_perturbation * sin(9.0 * magnetic_kx * x);  
                magnetic_field_z(ix) +=     B_perturbation * sin(10.0 * magnetic_kx * x);  
                magnetic_field_z(ix) +=     B_perturbation * sin(11.0 * magnetic_kx * x);  
                magnetic_field_z(ix) +=     B_perturbation * sin(12.0 * magnetic_kx * x);  
                magnetic_field_z(ix) +=     B_perturbation * sin(13.0 * magnetic_kx * x);  
                magnetic_field_z(ix) +=     B_perturbation * sin(14.0 * magnetic_kx * x);  
                magnetic_field_z(ix) +=     B_perturbation * sin(15.0 * magnetic_kx * x);  
                magnetic_field_z(ix) +=     B_perturbation * sin(16.0 * magnetic_kx * x);  
                
                magnetic_field_z(ix) +=     B_perturbation * sin(17.0 * magnetic_kx * x);  
                magnetic_field_z(ix) +=     B_perturbation * sin(18.0 * magnetic_kx * x);  
                magnetic_field_z(ix) +=     B_perturbation * sin(19.0 * magnetic_kx * x);  
                magnetic_field_z(ix) +=     B_perturbation * sin(20.0 * magnetic_kx * x);  
                magnetic_field_z(ix) +=     B_perturbation * sin(21.0 * magnetic_kx * x);  
                magnetic_field_z(ix) +=     B_perturbation * sin(22.0 * magnetic_kx * x);  
                magnetic_field_z(ix) +=     B_perturbation * sin(23.0 * magnetic_kx * x);  
                magnetic_field_z(ix) +=     B_perturbation * sin(24.0 * magnetic_kx * x);  
                magnetic_field_z(ix) +=     B_perturbation * sin(25.0 * magnetic_kx * x);  
                magnetic_field_z(ix) +=     B_perturbation * sin(26.0 * magnetic_kx * x);  
                magnetic_field_z(ix) +=     B_perturbation * sin(27.0 * magnetic_kx * x);  
                magnetic_field_z(ix) +=     B_perturbation * sin(28.0 * magnetic_kx * x);  
                magnetic_field_z(ix) +=     B_perturbation * sin(29.0 * magnetic_kx * x);  
                magnetic_field_z(ix) +=     B_perturbation * sin(30.0 * magnetic_kx * x);  
                magnetic_field_z(ix) +=     B_perturbation * sin(31.0 * magnetic_kx * x);  
                magnetic_field_z(ix) +=     B_perturbation * sin(32.0 * magnetic_kx * x);
                */
                



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
