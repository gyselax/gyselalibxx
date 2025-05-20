// SPDX-License-Identifier: MIT

#include <cassert>
#include <cmath>
#include <complex>
#include <iostream>

#include <ddc/ddc.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "geometry.hpp"
#include "hybridfieldsolver.hpp"
#include "vector_index_tools.hpp"

HybridFieldSolver::HybridFieldSolver(HybridSolver const& solve_hybrid, IMomentsCalculator const& compute_moments)
    : m_solve_hybrid(solve_hybrid)
    , m_compute_moments(compute_moments)
{
}

void HybridFieldSolver::operator()(
        DFieldXY const electrostatic_potential,
        DFieldXY const electric_field_x,
        DFieldXY const electric_field_y,
        DConstFieldSpVxVyXY const allfdistribu) const
{
    Kokkos::Profiling::pushRegion("HybridFieldSolver");
    assert((get_idx_range(electrostatic_potential) == get_idx_range<GridX, GridY>(allfdistribu)));
    IdxRangeXY const idx_range_xy = get_idx_range(electrostatic_potential);

    // Compute the RHS of the Quasi-Neutrality equation.
    DFieldMemXY rho(idx_range_xy);
    DFieldMemVxVy contiguous_slice_vxvy(get_idx_range<GridVx, GridVy>(allfdistribu));
    m_compute_moments(get_field(rho), allfdistribu);

    VectorField<
            double,
            IdxRangeXY,
            VectorIndexSet<X, Y>,
            Kokkos::DefaultExecutionSpace::memory_space,
            typename DFieldMemXY::layout_type>
            electric_field(electric_field_x, electric_field_y);
    m_solve_hybrid(electrostatic_potential, electric_field, get_field(rho));

    Kokkos::Profiling::popRegion();
}


void HybridFieldSolver::operator()(
    DFieldXY const pressure_field,
    DFieldXY const pressure_field_old,
    DFieldXY const pressure_field_mid,
    DFieldXY const pressure_field_previous,
    DFieldXY const magnetic_field_z,
    DFieldXY const magnetic_field_z_old,
    DFieldXY const magnetic_field_z_mid,
    DFieldXY const magnetic_field_z_previous,
    DFieldXY const mean_current_x_mid,
    DFieldXY const mean_current_y_mid,
    DFieldXY const rhs_1,
    DFieldXY const rhs_2,
    DFieldXY const rhs_3,
    DFieldXY const rhs_4,
    DFieldXY const mean_current_x,
    DFieldXY const mean_current_y,
    DFieldSpXY const mean_current_x_each,
    DFieldSpXY mean_current_y_each,
    DFieldSpXY rho_each,
    DFieldXY gradx_magnetic,
    DFieldXY grady_magnetic,
    DFieldXY gradx_pressure,
    DFieldXY grady_pressure,
    DFieldXY const rho,
    DConstFieldSpVxVyXY const allfdistribu,
    double const dt) const
{
Kokkos::Profiling::pushRegion("Hybridfieldsolver");
//assert((get_idx_range(electrostatic_potential) == get_idx_range<GridX, GridY>(allfdistribu)));

    m_compute_moments(rho, allfdistribu);
    m_compute_moments(mean_current_x, mean_current_y, rho, allfdistribu);

    m_compute_moments(rho_each, allfdistribu);
    m_compute_moments(mean_current_x_each, mean_current_y_each, rho_each, allfdistribu);

    IdxRangeXY idx_range(get_idx_range(magnetic_field_z));
    /*
    IdxRangeSp const kin_species_idx_range = get_idx_range<Species>(allfdistribu);
    
        for (IdxSp isp : kin_species_idx_range) {
            ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                idx_range,
                KOKKOS_LAMBDA(IdxXY const ixy) {
                    std::cout << "rho is: " << std::setprecision(16) << mean_velocity_y(ixy) << std::endl;
                    std::cout << "rho each is: " << std::setprecision(16) << mean_velocity_y_each(isp, ixy) << std::endl;
                });
        }
    */
    
            

    m_solve_hybrid(pressure_field, pressure_field_old, pressure_field_mid, pressure_field_previous, 
                magnetic_field_z, magnetic_field_z_old, magnetic_field_z_mid, magnetic_field_z_previous, 
                mean_current_x_mid, mean_current_y_mid, rhs_1, rhs_2, rhs_3, rhs_4,
                mean_current_x, mean_current_y, mean_current_x_each, mean_current_y_each, 
                rho_each, gradx_magnetic, grady_magnetic, gradx_pressure, grady_pressure, rho, dt);

Kokkos::Profiling::popRegion();
}
