// SPDX-License-Identifier: MIT

#include <cassert>
#include <cmath>
#include <complex>
#include <iostream>

#include <ddc/ddc.hpp>

#include <geometry.hpp>

#include "qnsolver.hpp"

QNSolver::QNSolver(PoissonSolver const& solve_poisson, IChargeDensityCalculator const& compute_rho)
    : m_solve_poisson(solve_poisson)
    , m_compute_rho(compute_rho)
{
}

void QNSolver::operator()(
        DSpanX const electrostatic_potential,
        DSpanX const electric_field,
        DViewSpXVx const allfdistribu) const
{
    Kokkos::Profiling::pushRegion("QNSolver");
    assert(electrostatic_potential.domain() == ddc::get_domain<IDimX>(allfdistribu));
    IDomainX const x_dom = electrostatic_potential.domain();
    // Compute the RHS of the Quasi-Neutrality equation.
    DFieldX rho(x_dom);

    m_compute_rho(rho, allfdistribu);

    m_solve_poisson(electrostatic_potential, electric_field, rho);

    Kokkos::Profiling::popRegion();
}
