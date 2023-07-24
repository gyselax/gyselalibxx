#include <stdexcept>
#include <string>

#include <ddc/ddc.hpp>

#include <maxwellianequilibrium.hpp>

#include "krook_source_constant.hpp"
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
 * amplitude is a constant
 *
 * therefore : 
 * f(t+dt) = ftarget + (f(t)-ftarget)*exp(-amplitude*mask*dt)
 */
KrookSourceConstant::KrookSourceConstant(
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
        ddc::expose_to_pdi("krook_source_constant_extent", m_extent);
        ddc::expose_to_pdi("krook_source_constant_stiffness", m_stiffness);
        ddc::expose_to_pdi("krook_source_constant_amplitude", m_amplitude);
        ddc::expose_to_pdi("krook_source_constant_density", m_density);
        ddc::expose_to_pdi("krook_source_constant_temperature", m_temperature);
        ddc::expose_to_pdi("krook_source_constant_ftarget", m_ftarget);
        ddc::expose_to_pdi("krook_source_constant_mask", m_mask);
        break;
    case RhsType::Sink:
        ddc::expose_to_pdi("krook_sink_constant_extent", m_extent);
        ddc::expose_to_pdi("krook_sink_constant_stiffness", m_stiffness);
        ddc::expose_to_pdi("krook_sink_constant_amplitude", m_amplitude);
        ddc::expose_to_pdi("krook_sink_constant_density", m_density);
        ddc::expose_to_pdi("krook_sink_constant_temperature", m_temperature);
        ddc::expose_to_pdi("krook_sink_constant_ftarget", m_ftarget);
        ddc::expose_to_pdi("krook_sink_constant_mask", m_mask);
        break;
    }
}

DSpanSpXVx KrookSourceConstant::operator()(DSpanSpXVx const allfdistribu, double const dt) const
{
    ddc::for_each(
            ddc::policies::parallel_host,
            allfdistribu.domain(),
            [&](IndexSpXVx const ispxvx) {
                allfdistribu(ispxvx)
                        = m_ftarget(ddc::select<IDimVx>(ispxvx))
                          + (allfdistribu(ispxvx) - m_ftarget(ddc::select<IDimVx>(ispxvx)))
                                    * std::exp(
                                            -m_amplitude * m_mask(ddc::select<IDimX>(ispxvx)) * dt);
            });

    return allfdistribu;
}
