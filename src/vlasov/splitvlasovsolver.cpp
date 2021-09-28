#include "iadvectionvx.h"
#include "iadvectionx.h"
#include "splitvlasovsolver.h"

SplitVlasovSolver::SplitVlasovSolver(const IAdvectionX& advec_x, const IAdvectionVx& m_advec_vx)
    : m_advec_x(advec_x)
    , m_advec_vx(m_advec_vx)
{
}


DSpanXVx SplitVlasovSolver::operator()(
        DSpanXVx fdistribu,
        DViewX efield,
        double sqrt_me_on_mspecies,
        double dt) const
{
    m_advec_x(fdistribu, sqrt_me_on_mspecies, dt / 2);
    m_advec_vx(fdistribu, efield, sqrt_me_on_mspecies, dt);
    m_advec_x(fdistribu, sqrt_me_on_mspecies, dt / 2);
    return fdistribu;
}
