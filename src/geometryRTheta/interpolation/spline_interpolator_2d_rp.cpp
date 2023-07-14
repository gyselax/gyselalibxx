
#include <sll/spline_builder_2d.hpp>
#include <sll/spline_evaluator_2d.hpp>

#include "i_interpolator_2d_rp.hpp"
#include "spline_interpolator_2d_rp.hpp"


DSpanRP SplineInterpolatorRP::operator()(
        DSpanRP const inout_data,
        ddc::ChunkSpan<CoordRP const, IDomainRP> const coordinates) const
{
#ifndef NDEBUG
    // To ensure that the interpolator is C0, we ensure that
    // the value at (r=0,theta) is the same for all theta.
    auto r_domain = ddc::get_domain<IDimR>(inout_data);
    auto theta_domain = ddc::get_domain<IDimP>(inout_data);
    if (ddc::coordinate(r_domain.front()) == 0) {
        ddc::for_each(theta_domain, [&](IndexP const ip) {
            assert(("Unicity of the value at the center point:",
                    inout_data(r_domain.front(), ip)
                            == inout_data(r_domain.front(), theta_domain.front())));
        });
    }
#endif

    m_builder(m_coefs, inout_data);
    m_evaluator(inout_data, coordinates, m_coefs);
    return inout_data;
}
