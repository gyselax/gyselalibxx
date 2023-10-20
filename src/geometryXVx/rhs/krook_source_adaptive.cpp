#include <cassert>
#include <stdexcept>
#include <string>

#include <ddc/ddc.hpp>

#include <maxwellianequilibrium.hpp>
#include <quadrature.hpp>
#include <rk2.hpp>
#include <species_info.hpp>
#include <trapezoid_quadrature.hpp>

#include "krook_source_adaptive.hpp"
#include "mask_tanh.hpp"

/**
 * Solves the equation \f$\partial f / \partial_t\f$ = -amplitude * mask * ( f - ftarget) (BGK operator)
 * 
 * mask defines the spatial region where the operator is active. 
 * If type = Source, the mask equals one in the central zone of the plasma of width extent; 
 * If type = Sink, the mask equals zero in the central zone of the plasma of width extent;
 *
 * ftarget is a maxwellian characterized by density and temperature, and a zero fluid velocity.
 *
 * amplitude depends on space, time and the considered species so that: 
 * amplitude(ions) = m_amplitude = constant
 * amplitude(electrons, x, t) = m_amplitude
 *                  * (density_ions(x,t) - m_density) / (density_electrons(x,t) - m_density)
 * so that the operator conserves locally the charge. 
 */
KrookSourceAdaptive::KrookSourceAdaptive(
        IDomainX const& gridx,
        IDomainVx const& gridvx,
        RhsType const type,
        double const extent,
        double const stiffness,
        double const amplitude,
        double const density,
        double const temperature)
    : m_type(type)
    , m_extent(extent)
    , m_stiffness(stiffness)
    , m_amplitude(amplitude)
    , m_density(density)
    , m_temperature(temperature)
    , m_ftarget(gridvx)
{
    // mask that defines the region where the operator is active
    switch (m_type) {
    case RhsType::Source:
        // the mask equals one in the interval [x_left, x_right]
        m_mask = mask_tanh(gridx, m_extent, m_stiffness, MaskType::Normal, false);
        break;
    case RhsType::Sink:
        // the mask equals zero in the center of the plasma
        m_mask = mask_tanh(gridx, m_extent, m_stiffness, MaskType::Inverted, false);
        break;
    }

    // target distribution function
    MaxwellianEquilibrium::compute_maxwellian(m_ftarget, m_density, m_temperature, 0.);

    switch (m_type) {
    case RhsType::Source:
        ddc::expose_to_pdi("krook_source_adaptive_extent", m_extent);
        ddc::expose_to_pdi("krook_source_adaptive_stiffness", m_stiffness);
        ddc::expose_to_pdi("krook_source_adaptive_amplitude", m_amplitude);
        ddc::expose_to_pdi("krook_source_adaptive_density", m_density);
        ddc::expose_to_pdi("krook_source_adaptive_temperature", m_temperature);
        ddc::expose_to_pdi("krook_source_adaptive_ftarget", m_ftarget);
        ddc::expose_to_pdi("krook_source_adaptive_mask", m_mask);
        break;
    case RhsType::Sink:
        ddc::expose_to_pdi("krook_sink_adaptive_extent", m_extent);
        ddc::expose_to_pdi("krook_sink_adaptive_stiffness", m_stiffness);
        ddc::expose_to_pdi("krook_sink_adaptive_amplitude", m_amplitude);
        ddc::expose_to_pdi("krook_sink_adaptive_density", m_density);
        ddc::expose_to_pdi("krook_sink_adaptive_temperature", m_temperature);
        ddc::expose_to_pdi("krook_sink_adaptive_ftarget", m_ftarget);
        ddc::expose_to_pdi("krook_sink_adaptive_mask", m_mask);
        break;
    }
}

/**
 * Computes the amplitude parameter for each species: 
 * amplitudes(ions) = m_amplitude = constant
 * amplitudes(electrons, x, t) = m_amplitude
 *                  * (density_ions(x,t) - m_density) / (density_electrons(x,t) - m_density)
 */
void KrookSourceAdaptive::get_amplitudes(DSpanSp amplitudes, DViewSpVx const allfdistribu) const
{
    IDomainSp const dom_sp(ddc::get_domain<IDimSp>(allfdistribu));
    assert(dom_sp.size() == 2);
    assert(charge(dom_sp.front()) * charge(dom_sp.back()) < 0);
    std::optional<IndexSp> iion_opt;
    for (IndexSp const isp : dom_sp) {
        if (charge(isp) > 0) {
            iion_opt = isp;
        }
    }
    IndexSp iion(iion_opt.value());
    amplitudes(iion) = m_amplitude;

    IDomainVx const gridvx = allfdistribu.domain<IDimVx>();
    Quadrature<IDimVx> const integrate_v(trapezoid_quadrature_coefficients(gridvx));
    double const density_ion = integrate_v(allfdistribu[iion]);
    double const density_electron = integrate_v(allfdistribu[ielec()]);

    amplitudes(ielec()) = m_amplitude * (density_ion - m_density) / (density_electron - m_density);
}

void KrookSourceAdaptive::get_derivative(
        DSpanSpXVx df,
        DViewSpXVx allfdistribu,
        DViewSpXVx allfdistribu_start) const
{
    IDomainSpVx sp_vx_dom = ddc::get_domain<IDimSp, IDimVx>(allfdistribu);
    IDomainSp sp_dom = ddc::get_domain<IDimSp>(allfdistribu);
    IDomainX x_dom = ddc::get_domain<IDimX>(allfdistribu);

    DFieldSp amplitudes(sp_dom);
    ddc::for_each(x_dom, [&](IndexX const ix) {
        DFieldSpVx allfdistribu_slice(allfdistribu[ix]);
        get_amplitudes(amplitudes.span_view(), allfdistribu_slice);
        ddc::for_each(sp_vx_dom, [&](IndexSpVx const ispvx) {
            IndexSp isp(ddc::select<IDimSp>(ispvx));
            IndexVx ivx(ddc::select<IDimVx>(ispvx));
            IndexSpXVx ispxvx(isp, ix, ivx);
            df(ispxvx)
                    = -m_mask(ix) * amplitudes(isp) * (allfdistribu_start(ispxvx) - m_ftarget(ivx));
        });
    });
}

DSpanSpXVx KrookSourceAdaptive::operator()(DSpanSpXVx const allfdistribu, double const dt) const
{
    RK2<DFieldSpXVx> timestepper(allfdistribu.domain());
    timestepper.update(allfdistribu, dt, [&](DSpanSpXVx df, DViewSpXVx f) {
        get_derivative(df, f, allfdistribu);
    });

    return allfdistribu;
}
