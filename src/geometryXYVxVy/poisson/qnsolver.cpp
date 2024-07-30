// SPDX-License-Identifier: MIT

#include <cassert>
#include <cmath>
#include <complex>
#include <iostream>

#include <ddc/ddc.hpp>

#include <directional_tag.hpp>
#include <geometry.hpp>

#include "qnsolver.hpp"

QNSolver::QNSolver(PoissonSolver const& solve_poisson, IChargeDensityCalculator const& compute_rho)
    : m_solve_poisson(solve_poisson)
    , m_compute_rho(compute_rho)
{
}

void QNSolver::operator()(
        DSpanXY const electrostatic_potential,
        DSpanXY const electric_field_x,
        DSpanXY const electric_field_y,
        DViewSpXYVxVy const allfdistribu) const
{
    Kokkos::Profiling::pushRegion("QNSolver");
    assert((electrostatic_potential.domain() == ddc::get_domain<IDimX, IDimY>(allfdistribu)));
    IDomainXY const xy_dom = electrostatic_potential.domain();

    // Compute the RHS of the Quasi-Neutrality equation.
    DFieldXY rho(xy_dom);
    DFieldVxVy contiguous_slice_vxvy(allfdistribu.domain<IDimVx, IDimVy>());
    m_compute_rho(rho, allfdistribu);

    VectorFieldSpan<
            double,
            IDomainXY,
            NDTag<RDimX, RDimY>,
            typename DFieldXY::layout_type,
            Kokkos::DefaultExecutionSpace::memory_space>
            electric_field(electric_field_x, electric_field_y);
    m_solve_poisson(electrostatic_potential, electric_field, rho);

    Kokkos::Profiling::popRegion();
}
