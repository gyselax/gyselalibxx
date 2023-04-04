#include <cassert>
#include <stdexcept>
#include <string>

#include <ddc/ddc.hpp>

#include <maxwellianequilibrium.hpp>
#include <quadrature.hpp>
#include <trapezoid_quadrature.hpp>

#include "krook_source_adaptive.hpp"
#include "mask_tanh.hpp"

/**
 * Solves the equation \partial f / \partial_t = -amplitude * mask * ( f - ftarget) (BGK operator)
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
        RhsSolver const solver_name,
        double const extent,
        double const stiffness,
        double const amplitude,
        double const density,
        double const temperature)
    : m_type(type)
    , m_solver_name(solver_name)
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

    // solver for the module
    switch (m_solver_name) {
    case RhsSolver::Rk2:
        m_solver = std::make_unique<RK2_solver>(
                [this](DSpanVx rhs_val,
                       DViewSpXVx allfdistribu,
                       double time,
                       IndexSpX const& ispx) { rhs(rhs_val, allfdistribu, time, ispx); });
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
 * amplitude(ions) = m_amplitude = constant
 * amplitude(electrons, x, t) = m_amplitude
 *                  * (density_ions(x,t) - m_density) / (density_electrons(x,t) - m_density)
 */
double KrookSourceAdaptive::get_amplitudes(DViewSpXVx const allfdistribu, IndexSpX const& ispx)
        const
{
    if (charge(ddc::select<IDimSp>(ispx)) >= 0.) {
        return m_amplitude;
    }

    assert(ddc::get_domain<IDimSp>(allfdistribu).size() == 2);
    IndexSp isp_ion;
    if (charge(ddc::select<IDimSp>(ispx))
        != charge(ddc::get_domain<IDimSp>(allfdistribu).front())) {
        isp_ion = ddc::get_domain<IDimSp>(allfdistribu).front();
    } else {
        isp_ion = ddc::get_domain<IDimSp>(allfdistribu).back();
    }

    // compute densities for each species
    IDomainVx const gridvx = allfdistribu.domain<IDimVx>();
    Quadrature<IDimVx> const integrate_v(trapezoid_quadrature_coefficients(gridvx));
    double const density_ion = integrate_v(allfdistribu[isp_ion][ddc::select<IDimX>(ispx)]);
    double const density_electron = integrate_v(allfdistribu[ddc::select<IDimSp, IDimX>(ispx)]);

    double const amplitude
            = m_amplitude * (density_ion - m_density) / (density_electron - m_density);
    return amplitude;
}

/**
 * right hand side of the equation \partial f / \partial_t = -amplitude * mask * ( f - ftarget)
 */
void KrookSourceAdaptive::rhs(
        DSpanVx const rhs,
        DViewSpXVx const allfdistribu,
        double const time,
        IndexSpX const& ispx) const
{
    double const amplitude = get_amplitudes(allfdistribu, ispx);
    ddc::for_each(
            ddc::policies::parallel_host,
            ddc::get_domain<IDimVx>(allfdistribu),
            [&](IndexVx const ivx) {
                IndexSpXVx ispxvx(ddc::select<IDimSp>(ispx), ddc::select<IDimX>(ispx), ivx);
                rhs(ivx) = -m_mask(ddc::select<IDimX>(ispx)) * amplitude
                           * (allfdistribu(ispxvx) - m_ftarget(ivx));
            });
}

DSpanSpXVx KrookSourceAdaptive::operator()(DSpanSpXVx const allfdistribu, double const dt) const
{
    return (*m_solver)(allfdistribu, dt);
}
