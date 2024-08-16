// SPDX-License-Identifier: MIT

#include <cassert>
#include <cmath>
#include <complex>
#include <iostream>

#include <ddc/ddc.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "directional_tag.hpp"
#include "geometry.hpp"
#include "qnsolver.hpp"

QNSolver::QNSolver(PoissonSolver const& solve_poisson, IChargeDensityCalculator const& compute_rho)
    : m_solve_poisson(solve_poisson)
    , m_compute_rho(compute_rho)
{
}

void QNSolver::operator()(
        DFieldXY const electrostatic_potential,
        DFieldXY const electric_field_x,
        DFieldXY const electric_field_y,
        DConstFieldSpXYVxVy const allfdistribu) const
{
    Kokkos::Profiling::pushRegion("QNSolver");
    assert((get_idx_range(electrostatic_potential) == get_idx_range<GridX, GridY>(allfdistribu)));
    IdxRangeXY const idx_range_xy = get_idx_range(electrostatic_potential);

    // Compute the RHS of the Quasi-Neutrality equation.
    DFieldMemXY rho(idx_range_xy);
    DFieldMemVxVy contiguous_slice_vxvy(get_idx_range<GridVx, GridVy>(allfdistribu));
    m_compute_rho(rho, allfdistribu);

    VectorField<
            double,
            IdxRangeXY,
            NDTag<X, Y>,
            typename DFieldMemXY::layout_type,
            Kokkos::DefaultExecutionSpace::memory_space>
            electric_field(electric_field_x, electric_field_y);
    m_solve_poisson(electrostatic_potential, electric_field, rho);

    Kokkos::Profiling::popRegion();
}
