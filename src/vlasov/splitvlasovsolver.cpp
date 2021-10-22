#include "iadvectionvx.hpp"
#include "iadvectionx.hpp"
#include "splitvlasovsolver.hpp"

SplitVlasovSolver::SplitVlasovSolver(IAdvectionX const& advec_x, IAdvectionVx const& m_advec_vx)
    : m_advec_x(advec_x)
    , m_advec_vx(m_advec_vx)
{
}


DSpanSpXVx SplitVlasovSolver::operator()(
        DSpanSpXVx allfdistribu,
        DViewX electric_potential,
        double dt) const
{
    m_advec_x(allfdistribu, dt / 2);
    m_advec_vx(allfdistribu, electric_potential, dt);
    m_advec_x(allfdistribu, dt / 2);
    return allfdistribu;
}
