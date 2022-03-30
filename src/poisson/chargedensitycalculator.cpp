// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include <species_info.hpp>

#include "chargedensitycalculator.hpp"

ChargeDensityCalculator::ChargeDensityCalculator(
        SpeciesInformation const& species_info,
        SplineVxBuilder const& spline_vx_builder,
        SplineEvaluator<BSplinesVx> const& spline_vx_evaluator)
    : m_species_info(species_info)
    , m_spline_vx_builder(spline_vx_builder)
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
    Chunk<double, BSDomainVx> vx_spline_coef(m_spline_vx_builder.spline_domain());

    IndexSp last_kin_species = allfdistribu.domain<IDimSp>().back();
    IndexSp last_species = m_species_info.charge().domain().back();
    double chargedens_adiabspecies = 0.;
    if (last_kin_species != last_species) {
        chargedens_adiabspecies = double(m_species_info.charge()(last_species));
    }

    for_each(rho.domain(), [&](IndexX const ix) {
        rho(ix) = chargedens_adiabspecies;
        for_each(get_domain<IDimSp>(allfdistribu), [&](IndexSp const isp) {
            deepcopy(f_vx_slice, allfdistribu[isp][ix]);
            m_spline_vx_builder(
                    vx_spline_coef.span_view(),
                    f_vx_slice.span_cview(),
                    &m_derivs_vxmin,
                    &m_derivs_vxmax);
            rho(ix) += m_species_info.charge()(isp)
                       * m_spline_vx_evaluator.integrate(vx_spline_coef.span_cview());
        });
    });
}
