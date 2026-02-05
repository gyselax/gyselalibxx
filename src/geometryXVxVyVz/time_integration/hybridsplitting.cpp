// SPDX-License-Identifier: MIT

#include <cmath>
#include <iostream>

#include <ddc/ddc.hpp>
#include <ddc/pdi.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ihybridfieldsolver.hpp"
#include "ihybridvlasovsolver.hpp"
#include "imomentscalculator.hpp"
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


DFieldSpVxVyVzX Hybridsplitting::operator()(
        DFieldSpVxVyVzX const allfdistribu_v3D_split,
        DFieldSpX mean_current_x_each, 
        DFieldSpX mean_current_y_each,
        DFieldSpX mean_current_z_each,
        DFieldX mean_current_x,
        DFieldX mean_current_y,
        DFieldX mean_current_z,
        DFieldX momentum_x,
        DFieldX momentum_y,
        DFieldX momentum_z,
        DFieldX magnetic_field_x,
        DFieldX magnetic_field_y,
        DFieldX magnetic_field_z,
        DFieldX pressure,
        DFieldSpX rho_each, 
        DFieldX rho, 
        DFieldX kinetic, 
        double const dt,
        int const steps) const
{
    
    IdxRangeSpXVxVyVz idx_range_v3D_split_output_layout(get_idx_range(allfdistribu_v3D_split));
    DFieldMemSpXVxVyVz allfdistribu_v3D_split_output_layout(idx_range_v3D_split_output_layout);
    auto allfdistribu_host_alloc
            = ddc::create_mirror_view(get_field(allfdistribu_v3D_split_output_layout));
    host_t<DFieldSpXVxVyVz> allfdistribu_host = get_field(allfdistribu_host_alloc);

    // electrostatic potential and electric field (depending only on x)
    DFieldMemX electrostatic_potential(get_idx_range<GridX>(allfdistribu_v3D_split));
    DFieldMemX electric_field_x(get_idx_range<GridX>(allfdistribu_v3D_split));
    DFieldMemX electric_field_y(get_idx_range<GridX>(allfdistribu_v3D_split));
    host_t<DFieldMemX> magnetic_field_x_host(
        get_idx_range<GridX>(allfdistribu_v3D_split));
    host_t<DFieldMemX> magnetic_field_y_host(
        get_idx_range<GridX>(allfdistribu_v3D_split));
    host_t<DFieldMemX> magnetic_field_z_host(
        get_idx_range<GridX>(allfdistribu_v3D_split));
    host_t<DFieldMemSpVxVy> density_spvxvy_host(get_idx_range<Species,GridVx,GridVy>(allfdistribu_v3D_split));

    // a 2D memory block of the same size as fdistribu
    DFieldMemSpVxVyVzX allfdistribu_half_t(get_idx_range(allfdistribu_v3D_split));

    // define the intermidiate variables used in the iteration solver of the fields.
    DFieldMemX pressure_old(get_idx_range<GridX>(allfdistribu_v3D_split));
    DFieldMemX pressure_mid(get_idx_range<GridX>(allfdistribu_v3D_split));
    DFieldMemX pressure_previous(get_idx_range<GridX>(allfdistribu_v3D_split));

    DFieldMemX temp_iter_1(get_idx_range<GridX>(allfdistribu_v3D_split));
    DFieldMemX temp_iter_2(get_idx_range<GridX>(allfdistribu_v3D_split));
    DFieldMemX temp_iter_3(get_idx_range<GridX>(allfdistribu_v3D_split));
    DFieldMemX temp_iter_4(get_idx_range<GridX>(allfdistribu_v3D_split));

    // used for the field solver 
    DFieldMemX magnetic_field_y_old(get_idx_range<GridX>(allfdistribu_v3D_split));
    DFieldMemX magnetic_field_y_mid(get_idx_range<GridX>(allfdistribu_v3D_split));
    DFieldMemX magnetic_field_y_previous(get_idx_range<GridX>(allfdistribu_v3D_split));
    DFieldMemX magnetic_field_z_old(get_idx_range<GridX>(allfdistribu_v3D_split));
    DFieldMemX magnetic_field_z_mid(get_idx_range<GridX>(allfdistribu_v3D_split));
    DFieldMemX magnetic_field_z_previous(get_idx_range<GridX>(allfdistribu_v3D_split));
    DFieldMemX u_bar_x(get_idx_range<GridX>(allfdistribu_v3D_split));
    DFieldMemX u_bar_y(get_idx_range<GridX>(allfdistribu_v3D_split));
    DFieldMemX u_bar_z(get_idx_range<GridX>(allfdistribu_v3D_split));
    DFieldMemX gradx_rho(get_idx_range<GridX>(allfdistribu_v3D_split));
    DFieldMemX gradx_magnetic_field_y_mid(get_idx_range<GridX>(allfdistribu_v3D_split));
    DFieldMemX gradx_magnetic_field_z_mid(get_idx_range<GridX>(allfdistribu_v3D_split));
    DFieldMemX rhs_1(get_idx_range<GridX>(allfdistribu_v3D_split));
    DFieldMemX rhs_2(get_idx_range<GridX>(allfdistribu_v3D_split));
    DFieldMemX rhs_3(get_idx_range<GridX>(allfdistribu_v3D_split));
    DFieldMemX curl_rhs_2(get_idx_range<GridX>(allfdistribu_v3D_split));
    DFieldMemX curl_rhs_3(get_idx_range<GridX>(allfdistribu_v3D_split));
    DFieldMemX Mxx(get_idx_range<GridX>(allfdistribu_v3D_split));
    DFieldMemX Mxy(get_idx_range<GridX>(allfdistribu_v3D_split));
    DFieldMemX Mxz(get_idx_range<GridX>(allfdistribu_v3D_split));
    DFieldMemX Myx(get_idx_range<GridX>(allfdistribu_v3D_split));
    DFieldMemX Myy(get_idx_range<GridX>(allfdistribu_v3D_split));
    DFieldMemX Myz(get_idx_range<GridX>(allfdistribu_v3D_split));
    DFieldMemX Mzx(get_idx_range<GridX>(allfdistribu_v3D_split));
    DFieldMemX Mzy(get_idx_range<GridX>(allfdistribu_v3D_split));
    DFieldMemX Mzz(get_idx_range<GridX>(allfdistribu_v3D_split));
    DFieldMemX weighted_u_x(get_idx_range<GridX>(allfdistribu_v3D_split));
    DFieldMemX weighted_u_y(get_idx_range<GridX>(allfdistribu_v3D_split));
    DFieldMemX weighted_u_z(get_idx_range<GridX>(allfdistribu_v3D_split));
    DFieldMemX weighted_p_para_x(get_idx_range<GridX>(allfdistribu_v3D_split));
    DFieldMemX weighted_p_para_y(get_idx_range<GridX>(allfdistribu_v3D_split));
    DFieldMemX weighted_p_para_z(get_idx_range<GridX>(allfdistribu_v3D_split));
    DFieldMemX qx(get_idx_range<GridX>(allfdistribu_v3D_split));
    DFieldMemX qy(get_idx_range<GridX>(allfdistribu_v3D_split));
    DFieldMemX qz(get_idx_range<GridX>(allfdistribu_v3D_split));
    DFieldMemX p_parallel_x(get_idx_range<GridX>(allfdistribu_v3D_split));
    DFieldMemX p_parallel_y(get_idx_range<GridX>(allfdistribu_v3D_split));
    DFieldMemX p_parallel_z(get_idx_range<GridX>(allfdistribu_v3D_split));

    DFieldMemSpXVxVy density_spxvxvy(get_idx_range<Species,GridX,GridVx,GridVy>(allfdistribu_v3D_split));
    DFieldMemSpVxVy density_spvxvy(get_idx_range<Species,GridVx,GridVy>(allfdistribu_v3D_split));

    DFieldMemSpX yl_x(get_idx_range(mean_current_x_each));
    DFieldMemSpX yl_y(get_idx_range(mean_current_x_each));
    DFieldMemSpX yl_z(get_idx_range(mean_current_x_each));

    DFieldMemSpX y2_x(get_idx_range(mean_current_x_each));
    DFieldMemSpX y2_y(get_idx_range(mean_current_x_each));
    DFieldMemSpX y2_z(get_idx_range(mean_current_x_each));

    DFieldMemSpX y3_x(get_idx_range(mean_current_x_each));
    DFieldMemSpX y3_y(get_idx_range(mean_current_x_each));
    DFieldMemSpX y3_z(get_idx_range(mean_current_x_each));

    DFieldMemSpX yr_x(get_idx_range(mean_current_x_each));
    DFieldMemSpX yr_y(get_idx_range(mean_current_x_each));
    DFieldMemSpX yr_z(get_idx_range(mean_current_x_each));

    DFieldMemSpX multi_para_tem(get_idx_range(mean_current_x_each));
    DFieldMemSpX multi_perp_tem(get_idx_range(mean_current_x_each));
    DFieldMemX single_para_tem(get_idx_range<GridX>(allfdistribu_v3D_split));
    DFieldMemX single_perp_tem(get_idx_range<GridX>(allfdistribu_v3D_split));
               
    
    double old_magnetic_energy = 0.0;
    double old_magnetic_energy_x = 0.0;
    double old_magnetic_energy_y = 0.0;
    double old_magnetic_energy_z = 0.0;

    double old_kinetic_energy = 0.0;
    double old_pressure_energy = 0.0;
    double old_thermal_energy = 0.0;
    double old_energy = 0.0; 
    double old_mass = 0.0;
    double old_momentum_x = 0.0;
    double old_momentum_y = 0.0;
    double old_momentum_z = 0.0;
    double electron_temperature = 0.;
    
    int iter = 0;
    for (; iter < steps; ++iter) {

        double const iter_time = iter * dt;

        if (iter == 0) { // just for tesing 
                // Computation of the density, mean velocity, and kinetic energy density.
                m_moments_calculator(get_field(rho), get_const_field(allfdistribu_v3D_split));
                m_moments_calculator(get_field(momentum_x), get_field(momentum_y), get_field(momentum_z), 
                                    get_const_field(allfdistribu_v3D_split));
                m_moments_calculator(get_field(weighted_u_x), get_field(weighted_u_y), get_field(weighted_u_z), get_field(rho), get_const_field(allfdistribu_v3D_split));
                m_moments_calculator(get_field(kinetic), get_const_field(allfdistribu_v3D_split), 'k');

                old_magnetic_energy_x = magnetic_energy_int(Kokkos::DefaultExecutionSpace(), get_const_field(magnetic_field_x));
                old_magnetic_energy_y = magnetic_energy_int(Kokkos::DefaultExecutionSpace(), get_const_field(magnetic_field_y));
                old_magnetic_energy_z = magnetic_energy_int(Kokkos::DefaultExecutionSpace(), get_const_field(magnetic_field_z));
                old_magnetic_energy = old_magnetic_energy_x + old_magnetic_energy_y + old_magnetic_energy_z;
                old_thermal_energy = thermal_energy_int(Kokkos::DefaultExecutionSpace(), get_const_field(rho), electron_temperature);
                //old_kinetic_energy = magnetic_energy_int(Kokkos::DefaultExecutionSpace(), get_const_field(weighted_u_x));
                //old_kinetic_energy += magnetic_energy_int(Kokkos::DefaultExecutionSpace(), get_const_field(weighted_u_y));
                //old_kinetic_energy += magnetic_energy_int(Kokkos::DefaultExecutionSpace(), get_const_field(weighted_u_z));

                old_kinetic_energy = function_int(Kokkos::DefaultExecutionSpace(), get_const_field(kinetic));
                old_pressure_energy = function_int(Kokkos::DefaultExecutionSpace(), get_const_field(pressure));
                old_energy = old_kinetic_energy + old_magnetic_energy + old_thermal_energy;
                std::cout << "initial old_magnetic_x_energy is: " << std::setprecision(16) << old_magnetic_energy_x << std::endl;
                std::cout << "initial old_magnetic_y_energy is: " << std::setprecision(16) << old_magnetic_energy_y << std::endl;
                std::cout << "initial old_magnetic_z_energy is: " << std::setprecision(16) << old_magnetic_energy_z << std::endl;
                std::cout << "initial old_magnetic_energy is: " << std::setprecision(16) << old_magnetic_energy << std::endl;
                std::cout << "initial old_kinetic_energy is: " << std::setprecision(16) << old_kinetic_energy << std::endl;
                std::cout << "initial thermal_energy is: " << std::setprecision(16) << old_thermal_energy << std::endl;
                std::cout << "initial old_pressure_energy is: " << std::setprecision(16) << old_pressure_energy << std::endl;
                std::cout << "initial old_energy is: " << std::setprecision(16) << old_energy << std::endl;
                old_mass = function_int(Kokkos::DefaultExecutionSpace(), get_const_field(rho));
                old_momentum_x = function_int(Kokkos::DefaultExecutionSpace(), get_const_field(momentum_x));
                old_momentum_y = function_int(Kokkos::DefaultExecutionSpace(), get_const_field(momentum_y));
                old_momentum_z = function_int(Kokkos::DefaultExecutionSpace(), get_const_field(momentum_z));
                
        }
        
        // advect in x over dt / 2
        m_vlasov_solver(get_field(allfdistribu_v3D_split), dt / 2);
        
        
        // computation of the fields
        m_hybrid_solver(
                get_field(magnetic_field_y),
                get_field(magnetic_field_y_old),
                get_field(magnetic_field_y_mid),
                get_field(magnetic_field_y_previous),
                get_field(magnetic_field_z),
                get_field(magnetic_field_z_old),
                get_field(magnetic_field_z_mid),
                get_field(magnetic_field_z_previous),
                get_field(mean_current_x_each),
                get_field(mean_current_y_each),
                get_field(mean_current_z_each),
                get_field(u_bar_x),
                get_field(u_bar_y),
                get_field(u_bar_z),
                get_field(rho),
                get_field(rho_each),
                get_field(magnetic_field_x),
                get_field(gradx_rho),
                get_field(gradx_magnetic_field_y_mid),
                get_field(gradx_magnetic_field_z_mid),
                get_field(rhs_1),
                get_field(rhs_2),
                get_field(rhs_3),
                get_field(curl_rhs_2),
                get_field(curl_rhs_3),
                get_field(Mxx), get_field(Mxy), get_field(Mxz), 
                get_field(Myx), get_field(Myy), get_field(Myz), 
                get_field(Mzx), get_field(Mzy), get_field(Mzz), 
                get_field(weighted_u_x), get_field(weighted_u_y), get_field(weighted_u_z), 
                get_field(weighted_p_para_x), get_field(weighted_p_para_y), get_field(weighted_p_para_z), 
                get_field(qx), get_field(qy), get_field(qz), 
                get_field(p_parallel_x), get_field(p_parallel_y), get_field(p_parallel_z),
                get_const_field(allfdistribu_v3D_split),
                electron_temperature,
                dt);
        

        // advect in v over dt
        
        m_vlasov_solver(get_field(allfdistribu_v3D_split),  
                        get_field(rhs_1), get_field(rhs_2), get_field(rhs_3), 
                        get_field(p_parallel_x), get_field(p_parallel_y), get_field(p_parallel_z), 
                        get_field(magnetic_field_x), get_field(magnetic_field_y_mid), get_field(magnetic_field_z_mid),
                        dt);
        
        
        
        
        
        // advect in x over dt / 2
        m_vlasov_solver(get_field(allfdistribu_v3D_split), dt / 2);
        
       
        transpose_layout(
                Kokkos::DefaultExecutionSpace(),
                get_field(allfdistribu_v3D_split_output_layout),
                get_const_field(allfdistribu_v3D_split));
      
        // Computation of the density, mean velocity, and kinetic energy density.
        m_moments_calculator(get_field(rho), get_const_field(allfdistribu_v3D_split));
        m_moments_calculator(get_field(momentum_x), get_field(momentum_y), get_field(momentum_z), 
                            get_const_field(allfdistribu_v3D_split));
        m_moments_calculator(get_field(kinetic), get_const_field(allfdistribu_v3D_split), 'k');
        // compute the integral of the distribution about vz
        m_moments_calculator(get_field(density_spxvxvy), get_const_field(allfdistribu_v3D_split));

        // compute the temperature 
        m_moments_calculator(get_field(rho_each), get_const_field(allfdistribu_v3D_split));
        m_moments_calculator(get_field(mean_current_x_each), get_field(mean_current_y_each), 
                        get_field(mean_current_z_each), get_field(rho_each), get_const_field(allfdistribu_v3D_split));

        m_moments_calculator(get_field(multi_para_tem), get_field(multi_perp_tem), 
                        get_const_field(rho_each), get_const_field(magnetic_field_x),  get_const_field(magnetic_field_y),  get_const_field(magnetic_field_z),
                        get_const_field(mean_current_x_each),  get_const_field(mean_current_y_each),  get_const_field(mean_current_z_each),
                        get_const_field(allfdistribu_v3D_split));
          
        m_hybrid_solver(
                get_field(multi_para_tem),
                get_field(multi_perp_tem),
                get_field(single_para_tem),
                get_field(single_perp_tem));
        double para_temperature = function_int(Kokkos::DefaultExecutionSpace(), get_const_field(single_para_tem)) / ddc::discrete_space<BSplinesX>().ncells();
        double perp_temperature = function_int(Kokkos::DefaultExecutionSpace(), get_const_field(single_perp_tem)) / ddc::discrete_space<BSplinesX>().ncells();
        std::cout << "parallel temperature: " << para_temperature << std::endl;
        std::cout << "perpendicular temperature: " << perp_temperature << std::endl;

        // conservation properties calculation
        double magnetic_energy_x = magnetic_energy_int(Kokkos::DefaultExecutionSpace(), get_const_field(magnetic_field_x));
        std::cout << "magnetic field x energy is: " << std::setprecision(16) << magnetic_energy_x << std::endl;
        double magnetic_energy_y = magnetic_energy_int(Kokkos::DefaultExecutionSpace(), get_const_field(magnetic_field_y));
        std::cout << "magnetic field y energy is: " << std::setprecision(16) << magnetic_energy_y << std::endl;
        double magnetic_energy_z = magnetic_energy_int(Kokkos::DefaultExecutionSpace(), get_const_field(magnetic_field_z));
        std::cout << "magnetic field z energy is: " << std::setprecision(16) << magnetic_energy_z << std::endl;

        double thermal_energy = thermal_energy_int(Kokkos::DefaultExecutionSpace(), get_const_field(rho), electron_temperature);
        std::cout << "thermal energy is: " << std::setprecision(16) << thermal_energy << std::endl;

        double magnetic_energy = magnetic_energy_x + magnetic_energy_y + magnetic_energy_z;

        std::cout << "magnetic field energy is: " << std::setprecision(16) << magnetic_energy << std::endl;

        double kinetic_energy = function_int(Kokkos::DefaultExecutionSpace(), get_const_field(kinetic));

        //double kinetic_energy = magnetic_energy_int(Kokkos::DefaultExecutionSpace(), get_const_field(weighted_u_x));
        //kinetic_energy += magnetic_energy_int(Kokkos::DefaultExecutionSpace(), get_const_field(weighted_u_y));
        //kinetic_energy += magnetic_energy_int(Kokkos::DefaultExecutionSpace(), get_const_field(weighted_u_z));
        std::cout << "kinetic energy is: " << std::setprecision(16) << kinetic_energy << std::endl;

        double pressure_energy = function_int(Kokkos::DefaultExecutionSpace(), get_const_field(pressure));
        std::cout << "pressure energy is: " << std::setprecision(16) << pressure_energy << std::endl;

        std::cout << "total energy is: " << std::setprecision(16) << kinetic_energy + magnetic_energy + thermal_energy << std::endl;

        std::cout << "total energy error is: " << std::setprecision(16) << kinetic_energy + magnetic_energy + thermal_energy - old_energy << std::endl;
        //std::cout << "kinetic energy error is: " << std::setprecision(16) << kinetic_energy - old_kinetic_energy << std::endl;

        double mass = function_int(Kokkos::DefaultExecutionSpace(), get_const_field(rho));
        std::cout << "mass is: " << mass << std::endl;
        std::cout << "mass error is: " << mass - old_mass << std::endl;
        double new_momentum_x = function_int(Kokkos::DefaultExecutionSpace(), get_const_field(momentum_x));
        std::cout << "momentum x is: " << new_momentum_x << std::endl;
        std::cout << "momentum x error is: " << new_momentum_x - old_momentum_x << std::endl;
        double new_momentum_y = function_int(Kokkos::DefaultExecutionSpace(), get_const_field(momentum_y));
        std::cout << "momentum y is: " << new_momentum_y << std::endl;
        std::cout << "momentum y error is: " << new_momentum_y - old_momentum_y << std::endl;
        double new_momentum_z = function_int(Kokkos::DefaultExecutionSpace(), get_const_field(momentum_z));
        std::cout << "momentum z is: " << new_momentum_z << std::endl;
        std::cout << "momentum z error is: " << new_momentum_z - old_momentum_z << std::endl;

        // copies necessary to PDI
        ddc::parallel_deepcopy(
                allfdistribu_host,
                get_const_field(allfdistribu_v3D_split_output_layout));
        ddc::parallel_deepcopy(magnetic_field_x_host, magnetic_field_x);
        ddc::parallel_deepcopy(magnetic_field_y_host, magnetic_field_y);
        ddc::parallel_deepcopy(magnetic_field_z_host, magnetic_field_z);
        ddc::parallel_deepcopy(density_spvxvy_host, density_spvxvy);
        ddc::PdiEvent("iteration")
                .with("iter", iter)
                .with("time_saved", iter_time)
                .with("kinetic_energy", kinetic_energy)
                .with("magnetic_energy", magnetic_energy)
                .with("thermal_energy", thermal_energy)
                .with("parallel_temperature", para_temperature)
                .with("perpendicular_temperature", perp_temperature)
                .with("magnetic_field_x", magnetic_field_x_host)
                .with("magnetic_field_y", magnetic_field_y_host)
                .with("magnetic_field_z", magnetic_field_z_host)
                .with("densityspvxvy", density_spvxvy_host);

        std::cout << "time step just finished is: " << iter << std::endl;
    }
    

    double const final_time = iter * dt;

    transpose_layout(
            Kokkos::DefaultExecutionSpace(),
            get_field(allfdistribu_v3D_split_output_layout),
            get_const_field(allfdistribu_v3D_split));
    //copies necessary to PDI
    ddc::parallel_deepcopy(
            allfdistribu_host,
            get_const_field(allfdistribu_v3D_split_output_layout));
    
    ddc::parallel_deepcopy(magnetic_field_x_host, magnetic_field_x);
    ddc::parallel_deepcopy(magnetic_field_y_host, magnetic_field_y);
    ddc::parallel_deepcopy(magnetic_field_z_host, magnetic_field_z);
    ddc::parallel_deepcopy(density_spvxvy_host, density_spvxvy);
    ddc::PdiEvent("last_iteration")
            .with("iter", iter)
            .with("time_saved", final_time)
            .with("magnetic_field_x", magnetic_field_x_host)
            .with("magnetic_field_y", magnetic_field_y_host)
            .with("magnetic_field_z", magnetic_field_z_host)
            .with("densityspvxvy", density_spvxvy_host);
    //std::cout << "just check-----------------: " << std::setprecision(16) << final_time << std::endl;

    return allfdistribu_v3D_split;

}

