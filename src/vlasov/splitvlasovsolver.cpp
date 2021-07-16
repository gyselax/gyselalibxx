#include "iadvectionvx.h"
#include "iadvectionx.h"
#include "splitvlasovsolver.h"

SplitVlasovSolver::SplitVlasovSolver(const IAdvectionX& advec_x, const IAdvectionVx& m_advec_vx)
    : m_advec_x(advec_x)
    , m_advec_vx(m_advec_vx)
{
}


DBlockSpanXVx SplitVlasovSolver::operator()(DBlockSpanXVx fdistribu, double mass_ratio, double dt)
        const
{
    m_advec_x(fdistribu, mass_ratio, dt / 2);
    m_advec_vx(fdistribu, mass_ratio, dt);
    m_advec_x(fdistribu, mass_ratio, dt / 2);
    return fdistribu;
}
