// SPDX-License-Identifier: MIT

#include <cassert>
#include <cmath>
#include <complex>
#include <iostream>

#include <ddc/ddc.hpp>

#include "geometry.hpp"
#include "qnsolver.hpp"

QNSolver::QNSolver(PoissonSolver const& solve_poisson, IChargeDensityCalculator const& compute_rho)
    : m_solve_poisson(solve_poisson)
    , m_compute_rho(compute_rho)
{
}

void QNSolver::operator()(
        DFieldX const electrostatic_potential,
        DFieldX const electric_field,
        DConstFieldSpXVx const allfdistribu) const
{
    Kokkos::Profiling::pushRegion("QNSolver");
    assert(get_idx_range(electrostatic_potential) == get_idx_range<GridX>(allfdistribu));
    IdxRangeX const idx_range_x = get_idx_range(electrostatic_potential);
    // Compute the RHS of the Quasi-Neutrality equation.
    DFieldMemX rho(idx_range_x);

    m_compute_rho(get_field(rho), allfdistribu);

    m_solve_poisson(electrostatic_potential, electric_field, get_field(rho));

    Kokkos::Profiling::popRegion();
}
