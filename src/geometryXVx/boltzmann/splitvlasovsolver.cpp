// SPDX-License-Identifier: MIT

#include "iadvectionvx.hpp"
#include "iadvectionx.hpp"
#include "splitvlasovsolver.hpp"

SplitVlasovSolver::SplitVlasovSolver(
        IAdvectionSpatial<GeometryXVx, IDimX> const& advec_x,
        IAdvectionVelocity<GeometryXVx, IDimVx> const& advec_vx)
    : m_advec_x(advec_x)
    , m_advec_vx(advec_vx)
{
}


device_t<DSpanSpXVx> SplitVlasovSolver::operator()(
        device_t<DSpanSpXVx> const allfdistribu_device,
        DViewX const electric_field,
        double const dt) const
{
    auto allfdistribu_alloc = ddc::create_mirror_view_and_copy(allfdistribu_device);
    ddc::ChunkSpan allfdistribu = allfdistribu_alloc.span_view();
    m_advec_x(allfdistribu, dt / 2);
    m_advec_vx(allfdistribu, electric_field, dt);
    m_advec_x(allfdistribu, dt / 2);
    ddc::deepcopy(allfdistribu_device, allfdistribu);
    return allfdistribu_device;
}
