#include "iadvectionvx.h"
#include "iadvectionx.h"
#include "splitvlasovsolver.h"

SplitVlasovSolver::SplitVlasovSolver(const IAdvectionX& advec_x, const IAdvectionVx& m_advec_vx)
    : m_advec_x(advec_x)
    , m_advec_vx(m_advec_vx)
{
}


DSpanSpXVx SplitVlasovSolver::operator()(DSpanSpXVx fdistribu, DViewX electric_potential, double dt)
        const
{
    m_advec_x(fdistribu, dt / 2);
    m_advec_vx(fdistribu, electric_potential, dt);
    m_advec_x(fdistribu, dt / 2);
    return fdistribu;
}
