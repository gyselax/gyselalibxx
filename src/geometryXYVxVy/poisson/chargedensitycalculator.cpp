// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include "chargedensitycalculator.hpp"

ChargeDensityCalculator::ChargeDensityCalculator(
        SplineVxVyBuilder const& spline_vxvy_builder,
        SplineVxVyEvaluator const& spline_vxvy_evaluator)
    : m_spline_vxvy_builder(spline_vxvy_builder)
    , m_spline_vxvy_evaluator(spline_vxvy_evaluator)
    , m_nbc_Vx(BSplinesVx::degree() / 2)
    , m_nbc_Vy(BSplinesVy::degree() / 2)
    , m_interp_dom_size_Vx(m_spline_vxvy_builder.interpolation_domain1().size())
    , m_interp_dom_size_Vy(m_spline_vxvy_builder.interpolation_domain2().size())
{
}

void ChargeDensityCalculator::operator()(DSpanXY const rho, DViewSpXYVxVy const allfdistribu) const
{
    DFieldVxVy f_vxvy_slice(allfdistribu.domain<IDimVx, IDimVy>());
    ddc::Chunk<double, BSDomainVxVy> vxvy_spline_coef(m_spline_vxvy_builder.spline_domain());

    IndexSp const last_kin_species = allfdistribu.domain<IDimSp>().back();
    IndexSp const last_species = ddc::discrete_space<IDimSp>().charges().domain().back();
    double chargedens_adiabspecies = 0.;
    if (last_kin_species != last_species) {
        chargedens_adiabspecies = double(charge(last_species));
    }

    // Compute the derivatives at boundaries
    // TODO: Compute the derivatives instead of imposing value = 0 (see 2d_spline_builder.cpp)
    std::vector<double> deriv_vxmin_data_Vx(m_nbc_Vx * m_interp_dom_size_Vy, 0.);
    std::vector<double> deriv_vxmax_data_Vx(m_nbc_Vx * m_interp_dom_size_Vy, 0.);
    DSpan2D deriv_vxmin(deriv_vxmin_data_Vx.data(), m_interp_dom_size_Vy, m_nbc_Vx);
    DSpan2D deriv_vxmax(deriv_vxmax_data_Vx.data(), m_interp_dom_size_Vy, m_nbc_Vx);

    std::vector<double> deriv_vymin_data_Vy(m_nbc_Vy * m_interp_dom_size_Vx, 0.);
    std::vector<double> deriv_vymax_data_Vy(m_nbc_Vy * m_interp_dom_size_Vx, 0.);
    DSpan2D deriv_vymin(deriv_vymin_data_Vy.data(), m_interp_dom_size_Vx, m_nbc_Vy);
    DSpan2D deriv_vymax(deriv_vymax_data_Vy.data(), m_interp_dom_size_Vx, m_nbc_Vy);

    std::vector<double> mixed_deriv_vxmin_vymin_data(m_nbc_Vx * m_nbc_Vy);
    std::vector<double> mixed_deriv_vxmin_vymax_data(m_nbc_Vx * m_nbc_Vy);
    std::vector<double> mixed_deriv_vxmax_vymin_data(m_nbc_Vx * m_nbc_Vy);
    std::vector<double> mixed_deriv_vxmax_vymax_data(m_nbc_Vx * m_nbc_Vy);
    DSpan2D mixed_deriv_vxmin_vymin(mixed_deriv_vxmin_vymin_data.data(), m_nbc_Vx, m_nbc_Vy);
    DSpan2D mixed_deriv_vxmin_vymax(mixed_deriv_vxmin_vymax_data.data(), m_nbc_Vx, m_nbc_Vy);
    DSpan2D mixed_deriv_vxmax_vymin(mixed_deriv_vxmax_vymin_data.data(), m_nbc_Vx, m_nbc_Vy);
    DSpan2D mixed_deriv_vxmax_vymax(mixed_deriv_vxmax_vymax_data.data(), m_nbc_Vx, m_nbc_Vy);

    ddc::for_each(rho.domain(), [&](IndexXY const ixy) {
        IndexX const ix = ddc::select<IDimX>(ixy);
        IndexY const iy = ddc::select<IDimY>(ixy);
        rho(ix, iy) = chargedens_adiabspecies;
        ddc::for_each(ddc::get_domain<IDimSp>(allfdistribu), [&](IndexSp const isp) {
            ddc::deepcopy(f_vxvy_slice, allfdistribu[isp][ix][iy]);
            m_spline_vxvy_builder(
                    vxvy_spline_coef.span_view(),
                    f_vxvy_slice.span_cview(),
                    deriv_vxmin,
                    deriv_vxmax,
                    deriv_vymin,
                    deriv_vymax,
                    mixed_deriv_vxmin_vymin,
                    mixed_deriv_vxmin_vymax,
                    mixed_deriv_vxmax_vymin,
                    mixed_deriv_vxmax_vymax);
            rho(ix, iy) += charge(isp)
                           * m_spline_vxvy_evaluator.integrate(vxvy_spline_coef.span_cview());
        });
    });
}
