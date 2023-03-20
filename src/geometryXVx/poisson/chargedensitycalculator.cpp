// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include "chargedensitycalculator.hpp"

ChargeDensityCalculator::ChargeDensityCalculator(
        SplineVxBuilder const& spline_vx_builder,
        SplineEvaluator<BSplinesVx> const& spline_vx_evaluator)
    : m_spline_vx_builder(spline_vx_builder)
    , m_spline_vx_evaluator(spline_vx_evaluator)
    , m_derivs_vxmin_data(BSplinesVx::degree() / 2, 0.)
    , m_derivs_vxmin(m_derivs_vxmin_data.data(), m_derivs_vxmin_data.size())
    , m_derivs_vxmax_data(BSplinesVx::degree() / 2, 0.)
    , m_derivs_vxmax(m_derivs_vxmax_data.data(), m_derivs_vxmax_data.size())
{
}

void ChargeDensityCalculator::operator()(DSpanX const rho, DViewSpXVx const allfdistribu) const
{
    DFieldVx f_vx_slice(allfdistribu.domain<IDimVx>());
    ddc::Chunk<double, BSDomainVx> vx_spline_coef(m_spline_vx_builder.spline_domain());

    IndexSp const last_kin_species = allfdistribu.domain<IDimSp>().back();
    IndexSp const last_species = ddc::discrete_space<IDimSp>().charges().domain().back();
    double chargedens_adiabspecies = 0.;
    if (last_kin_species != last_species) {
        chargedens_adiabspecies = double(charge(last_species));
    }

    ddc::for_each(rho.domain(), [&](IndexX const ix) {
        rho(ix) = chargedens_adiabspecies;
        ddc::for_each(ddc::get_domain<IDimSp>(allfdistribu), [&](IndexSp const isp) {
            ddc::deepcopy(f_vx_slice, allfdistribu[isp][ix]);
            m_spline_vx_builder(
                    vx_spline_coef.span_view(),
                    f_vx_slice.span_cview(),
                    m_derivs_vxmin,
                    m_derivs_vxmax);
            rho(ix) += charge(isp) * m_spline_vx_evaluator.integrate(vx_spline_coef.span_cview());
        });
    });
}
