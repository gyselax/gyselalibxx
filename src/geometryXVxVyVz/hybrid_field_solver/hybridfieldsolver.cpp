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
    DFieldX magnetic_field_y, 
    DFieldX magnetic_field_y_old,               
    DFieldX magnetic_field_y_mid, 
    DFieldX magnetic_field_y_previous,
    DFieldX magnetic_field_z, 
    DFieldX magnetic_field_z_old, 
    DFieldX magnetic_field_z_mid, 
    DFieldX magnetic_field_z_previous,
    DFieldSpX u_old_x, 
    DFieldSpX u_old_y, 
    DFieldSpX u_old_z, 
    DFieldX u_bar_x, 
    DFieldX u_bar_y, 
    DFieldX u_bar_z, 
    DFieldX rho,
    DFieldSpX rho_each,
    DFieldX magnetic_field_x,
    DFieldX gradx_rho,
    DFieldX gradx_magnetic_field_y_mid,
    DFieldX gradx_magnetic_field_z_mid,
    DFieldX rhs_1,
    DFieldX rhs_2,
    DFieldX rhs_3,
    DFieldX rhs_5,
    DFieldX rhs_6,
    DFieldX Mxx, DFieldX Mxy, DFieldX Mxz,
    DFieldX Myx, DFieldX Myy, DFieldX Myz,
    DFieldX Mzx, DFieldX Mzy, DFieldX Mzz,
    DFieldX weighted_u_x, DFieldX weighted_u_y, DFieldX weighted_u_z,
    DFieldX weighted_p_para_x, DFieldX weighted_p_para_y, DFieldX weighted_p_para_z,
    DFieldX qx, DFieldX qy, DFieldX qz,
    DFieldX p_parallel_x, DFieldX p_parallel_y, DFieldX p_parallel_z,
    DConstFieldSpVxVyVzX const allfdistribu,
    double const electron_temperature,
    double const dt) const
{
Kokkos::Profiling::pushRegion("Hybridfieldsolver1d3v");
//assert((get_idx_range(electrostatic_potential) == get_idx_range<GridX, GridY>(allfdistribu)));

    m_compute_moments(rho, allfdistribu);
    //m_compute_moments(weighted_u_x, weighted_u_y, weighted_u_z, rho, allfdistribu);

    m_compute_moments(rho_each, allfdistribu);
    m_compute_moments(u_old_x, u_old_y, u_old_z, rho_each, allfdistribu);

    IdxRangeX idx_range(get_idx_range(magnetic_field_z));
    
    m_solve_hybrid(
        magnetic_field_y, magnetic_field_y_old, magnetic_field_y_mid, magnetic_field_y_previous,
        magnetic_field_z, magnetic_field_z_old, magnetic_field_z_mid, magnetic_field_z_previous,
        u_old_x, u_old_y, u_old_z, 
        u_bar_x, u_bar_y, u_bar_z, 
        rho, rho_each,
        magnetic_field_x, gradx_rho, gradx_magnetic_field_y_mid, gradx_magnetic_field_z_mid,
        rhs_1, rhs_2, rhs_3, rhs_5, rhs_6,         
        Mxx, Mxy, Mxz, Myx, Myy, Myz, Mzx, Mzy, Mzz,
        weighted_u_x, weighted_u_y, weighted_u_z,
        weighted_p_para_x, weighted_p_para_y, weighted_p_para_z,
        qx, qy, qz, p_parallel_x, p_parallel_y, p_parallel_z,
        electron_temperature, dt      
        );    
    
Kokkos::Profiling::popRegion();
}



void HybridFieldSolver::operator()(
    DFieldSpX multi_para_tem, 
    DFieldSpX multi_perp_tem, 
    DFieldX single_para_tem, 
    DFieldX single_perp_tem
    ) const
{
Kokkos::Profiling::pushRegion("Hybridfieldsolver1d3v");
    m_solve_hybrid(
        multi_para_tem, multi_perp_tem, single_para_tem, single_perp_tem);    
    
Kokkos::Profiling::popRegion();
}

