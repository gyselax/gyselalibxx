#include "splitvlasovsolver.h"

SplitVlasovSolver::SplitVlasovSolver(const IAdvectionX& advec_x, const IAdvectionX& m_advec_vx)
    : m_advec_x(advec_x)
    , m_advec_vx(m_advec_vx)
{
}


DBlockViewXVx& SplitVlasovSolver::operator()(DBlockViewXVx& fdistribu, double mass_ratio, double dt) const
{
    m_advec_x(fdistribu, mass_ratio, dt / 2);
    m_advec_vx(fdistribu, mass_ratio, dt);
    m_advec_x(fdistribu, mass_ratio, dt / 2);
    return fdistribu;
}
