// SPDX-License-Identifier: MIT

#include <cmath>
#include <iostream>

#include <ddc/ddc.hpp>
#include <ddc/pdi.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ihybridfieldsolver.hpp"
#include "ihybridvlasovsolver.hpp"
#include "imomentscalculator.hpp"
#include "iqnsolver.hpp"
#include "hybridsplitting.hpp"
#include "transpose.hpp"
#include "l_norm_tools.hpp"

Hybridsplitting::Hybridsplitting(IHybridVlasovSolver const& vlasov_solver, IHybridFieldSolver const& hybrid_solver, 
                                IMomentsCalculator const& moments_calculator)
    : m_vlasov_solver(vlasov_solver)
    , m_hybrid_solver(hybrid_solver)
    , m_moments_calculator(moments_calculator)
{
}


DFieldSpVxVyXY Hybridsplitting::operator()(
        DFieldSpVxVyXY const allfdistribu_v2D_split,
        DFieldSpXY mean_current_x_each, 
        DFieldSpXY mean_current_y_each,
        DFieldXY mean_current_x,
        DFieldXY mean_current_y,
        DFieldXY momentum_x,
        DFieldXY momentum_y,
        DFieldXY magnetic_field_z,
        DFieldXY pressure,
        DFieldSpXY rho_each, 
        DFieldXY rho, 
        DFieldXY kinetic, 
        double const dt,
        int const steps) const
{
    
    IdxRangeSpXYVxVy idx_range_v2D_split_output_layout(get_idx_range(allfdistribu_v2D_split));
    DFieldMemSpXYVxVy allfdistribu_v2D_split_output_layout(idx_range_v2D_split_output_layout);
    auto allfdistribu_host_alloc
            = ddc::create_mirror_view(get_field(allfdistribu_v2D_split_output_layout));
    host_t<DFieldSpXYVxVy> allfdistribu_host = get_field(allfdistribu_host_alloc);

    // electrostatic potential and electric field (depending only on x)
    DFieldMemXY electrostatic_potential(get_idx_range<GridX, GridY>(allfdistribu_v2D_split));
    DFieldMemXY electric_field_x(get_idx_range<GridX, GridY>(allfdistribu_v2D_split));
    DFieldMemXY electric_field_y(get_idx_range<GridX, GridY>(allfdistribu_v2D_split));
    host_t<DFieldMemXY> electrostatic_potential_host(
        get_idx_range<GridX, GridY>(allfdistribu_v2D_split));

    // a 2D memory block of the same size as fdistribu
    DFieldMemSpVxVyXY allfdistribu_half_t(get_idx_range(allfdistribu_v2D_split));

    // define the intermidiate variables used in the iteration solver of the fields.
    DFieldMemXY gradx_magnetic(get_idx_range<GridX, GridY>(allfdistribu_v2D_split));
    DFieldMemXY grady_magnetic(get_idx_range<GridX, GridY>(allfdistribu_v2D_split));
    DFieldMemXY gradx_pressure(get_idx_range<GridX, GridY>(allfdistribu_v2D_split));
    DFieldMemXY grady_pressure(get_idx_range<GridX, GridY>(allfdistribu_v2D_split));

    DFieldMemXY pressure_old(get_idx_range<GridX, GridY>(allfdistribu_v2D_split));
    DFieldMemXY pressure_mid(get_idx_range<GridX, GridY>(allfdistribu_v2D_split));
    DFieldMemXY pressure_previous(get_idx_range<GridX, GridY>(allfdistribu_v2D_split));

    DFieldMemXY magnetic_field_old(get_idx_range<GridX, GridY>(allfdistribu_v2D_split));
    DFieldMemXY magnetic_field_mid(get_idx_range<GridX, GridY>(allfdistribu_v2D_split));
    DFieldMemXY magnetic_field_previous(get_idx_range<GridX, GridY>(allfdistribu_v2D_split));

    DFieldMemXY temp_iter_1(get_idx_range<GridX, GridY>(allfdistribu_v2D_split));
    DFieldMemXY temp_iter_2(get_idx_range<GridX, GridY>(allfdistribu_v2D_split));
    DFieldMemXY temp_iter_3(get_idx_range<GridX, GridY>(allfdistribu_v2D_split));
    DFieldMemXY temp_iter_4(get_idx_range<GridX, GridY>(allfdistribu_v2D_split));

    DFieldMemXY ux_bar(get_idx_range<GridX, GridY>(allfdistribu_v2D_split));
    DFieldMemXY uy_bar(get_idx_range<GridX, GridY>(allfdistribu_v2D_split));

    double old_magnetic_energy = 0.0;
    double old_kinetic_energy = 0.0;
    double old_pressure_energy = 0.0;
    double old_energy = 0.0; 
    double old_mass = 0.0;
    double old_momentum_x = 0.0;
    double old_momentum_y = 0.0;
    
    int iter = 0;
    for (; iter < steps; ++iter) {

        double const iter_time = iter * dt;

        if (iter == 0) { // just for tesing 
                // Computation of the density, mean velocity, and kinetic energy density.
                m_moments_calculator(get_field(rho), get_const_field(allfdistribu_v2D_split));
                m_moments_calculator(get_field(momentum_x), get_field(momentum_y), 
                                    get_const_field(allfdistribu_v2D_split));
                m_moments_calculator(get_field(kinetic), get_const_field(allfdistribu_v2D_split), 'k');

                old_magnetic_energy = magnetic_energy_int(Kokkos::DefaultExecutionSpace(), get_const_field(magnetic_field_z));
                old_kinetic_energy = function_int(Kokkos::DefaultExecutionSpace(), get_const_field(kinetic));
                old_pressure_energy = function_int(Kokkos::DefaultExecutionSpace(), get_const_field(pressure));
                old_energy = old_kinetic_energy + old_magnetic_energy + 1.5 * old_pressure_energy;
                std::cout << "initial old_magnetic_energy is: " << std::setprecision(16) << old_magnetic_energy << std::endl;
                std::cout << "initial old_kinetic_energy is: " << std::setprecision(16) << old_kinetic_energy << std::endl;
                std::cout << "initial old_pressure_energy is: " << std::setprecision(16) << old_pressure_energy << std::endl;
                std::cout << "initial old_energy is: " << std::setprecision(16) << old_energy << std::endl;
                old_mass = function_int(Kokkos::DefaultExecutionSpace(), get_const_field(rho));
                old_momentum_x = function_int(Kokkos::DefaultExecutionSpace(), get_const_field(momentum_x));
                old_momentum_y = function_int(Kokkos::DefaultExecutionSpace(), get_const_field(momentum_y));
                
        }

        // advect in x over dt / 2
        m_vlasov_solver(get_field(allfdistribu_v2D_split), dt / 2);
        //std::cout << "just check-----------------: " << std::setprecision(16) << iter_time << std::endl;

        // computation of the fields
        m_hybrid_solver(
                get_field(pressure),
                get_field(pressure_old),
                get_field(pressure_mid),
                get_field(pressure_previous),
                get_field(magnetic_field_z),
                get_field(magnetic_field_old),
                get_field(magnetic_field_mid),
                get_field(magnetic_field_previous),
                get_field(ux_bar),
                get_field(uy_bar),
                get_field(temp_iter_1),
                get_field(temp_iter_2),
                get_field(temp_iter_3),
                get_field(temp_iter_4),
                get_field(mean_current_x),
                get_field(mean_current_y),
                get_field(mean_current_x_each),
                get_field(mean_current_y_each),
                get_field(rho_each),
                get_field(gradx_magnetic),
                get_field(grady_magnetic),
                get_field(gradx_pressure),
                get_field(grady_pressure),
                get_field(rho),
                get_const_field(allfdistribu_v2D_split),
                dt);

        // advect in v over dt
        m_vlasov_solver(get_field(allfdistribu_v2D_split), get_field(magnetic_field_mid), 
                        get_field(temp_iter_1), get_field(temp_iter_2), dt);
        
        // advect in x over dt / 2
        m_vlasov_solver(get_field(allfdistribu_v2D_split), dt / 2);
        
        transpose_layout(
                Kokkos::DefaultExecutionSpace(),
                get_field(allfdistribu_v2D_split_output_layout),
                get_const_field(allfdistribu_v2D_split));
        // copies necessary to PDI
        ddc::parallel_deepcopy(
                allfdistribu_host,
                get_const_field(allfdistribu_v2D_split_output_layout));
        ddc::parallel_deepcopy(electrostatic_potential_host, magnetic_field_z);
        ddc::PdiEvent("iteration")
                .with("iter", iter)
                .with("time_saved", iter_time)
                .with("electrostatic_potential", electrostatic_potential_host);

        // Computation of the density, mean velocity, and kinetic energy density.
        m_moments_calculator(get_field(rho), get_const_field(allfdistribu_v2D_split));
        m_moments_calculator(get_field(momentum_x), get_field(momentum_y), 
                            get_const_field(allfdistribu_v2D_split));
        m_moments_calculator(get_field(kinetic), get_const_field(allfdistribu_v2D_split), 'k');

        std::cout << "time step just finished is: " << iter << std::endl;
        // conservation properties calculation
        double magnetic_energy = magnetic_energy_int(Kokkos::DefaultExecutionSpace(), get_const_field(magnetic_field_z));
        std::cout << "magnetic field energy is: " << std::setprecision(16) << magnetic_energy << std::endl;

        double kinetic_energy = function_int(Kokkos::DefaultExecutionSpace(), get_const_field(kinetic));
        std::cout << "kinetic energy is: " << std::setprecision(16) << kinetic_energy << std::endl;

        double pressure_energy = function_int(Kokkos::DefaultExecutionSpace(), get_const_field(pressure));
        std::cout << "pressure energy is: " << std::setprecision(16) << pressure_energy << std::endl;

        std::cout << "total energy error is: " << std::setprecision(16) << kinetic_energy + magnetic_energy + 1.5 * pressure_energy - old_energy << std::endl;
        //std::cout << "kinetic energy error is: " << std::setprecision(16) << kinetic_energy - old_kinetic_energy << std::endl;

        double mass = function_int(Kokkos::DefaultExecutionSpace(), get_const_field(rho));
        std::cout << "mass is: " << mass << std::endl;
        std::cout << "mass error is: " << mass - old_mass << std::endl;
        double new_momentum_x = function_int(Kokkos::DefaultExecutionSpace(), get_const_field(momentum_x));
        std::cout << "momentum x error is: " << new_momentum_x - old_momentum_x << std::endl;
        double new_momentum_y = function_int(Kokkos::DefaultExecutionSpace(), get_const_field(momentum_y));
        std::cout << "momentum y error is: " << new_momentum_y - old_momentum_y << std::endl;
    }
    

    double const final_time = iter * dt;

    m_hybrid_solver(
        get_field(electrostatic_potential),
        get_field(electric_field_x),
        get_field(electric_field_y),
        get_const_field(allfdistribu_v2D_split));

    transpose_layout(
            Kokkos::DefaultExecutionSpace(),
            get_field(allfdistribu_v2D_split_output_layout),
            get_const_field(allfdistribu_v2D_split));
    //copies necessary to PDI
    ddc::parallel_deepcopy(
            allfdistribu_host,
            get_const_field(allfdistribu_v2D_split_output_layout));
    
    ddc::parallel_deepcopy(electrostatic_potential_host, magnetic_field_z);
    ddc::PdiEvent("last_iteration")
            .with("iter", iter)
            .with("time_saved", final_time)
            .with("electrostatic_potential", electrostatic_potential_host);
    //std::cout << "just check-----------------: " << std::setprecision(16) << final_time << std::endl;

    return allfdistribu_v2D_split;

}

