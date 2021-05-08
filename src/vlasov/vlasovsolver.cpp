#include "vlasovsolver.h"

VlasovSolver::VlasovSolver(const Advection1D& advec_x, const Advection1D& m_advec_vx)
    : m_advec_x(advec_x)
    , m_advec_vx(m_advec_vx)
{
}


void VlasovSolver::operator()(DBlock2D& fdistribu, double mass_ratio, double dt) const
{
    m_advec_x(fdistribu, mass_ratio, dt / 2);
    m_advec_vx(fdistribu, mass_ratio, dt);
    m_advec_x(fdistribu, mass_ratio, dt / 2);
}
