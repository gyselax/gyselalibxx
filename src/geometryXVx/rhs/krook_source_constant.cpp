#include <stdexcept>
#include <string>

#include <ddc/ddc.hpp>

#include <maxwellianequilibrium.hpp>

#include "krook_source_constant.hpp"
#include "mask_tanh.hpp"

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
    , m_mask(gridx)
{
    // mask that defines the region where the operator is active
    DFieldX mask_host(gridx);
    switch (m_type) {
    case RhsType::Source:
        // the mask equals one in the interval [x_left, x_right]
        mask_host = mask_tanh(gridx, m_extent, m_stiffness, MaskType::Normal, false);
        break;
    case RhsType::Sink:
        // the mask equals zero in the center of the plasma
        mask_host = mask_tanh(gridx, m_extent, m_stiffness, MaskType::Inverted, false);
        break;
    }
    ddc::deepcopy(m_mask.span_view(), mask_host);

    // target distribution function
    MaxwellianEquilibrium::compute_maxwellian(m_ftarget.span_view(), m_density, m_temperature, 0.);
    auto ftarget_host = ddc::create_mirror_view_and_copy(m_ftarget.span_view());

    switch (m_type) {
    case RhsType::Source:
        ddc::expose_to_pdi("krook_source_constant_extent", m_extent);
        ddc::expose_to_pdi("krook_source_constant_stiffness", m_stiffness);
        ddc::expose_to_pdi("krook_source_constant_amplitude", m_amplitude);
        ddc::expose_to_pdi("krook_source_constant_density", m_density);
        ddc::expose_to_pdi("krook_source_constant_temperature", m_temperature);
        ddc::expose_to_pdi("krook_source_constant_ftarget", ftarget_host);
        ddc::expose_to_pdi("krook_source_constant_mask", mask_host);
        break;
    case RhsType::Sink:
        ddc::expose_to_pdi("krook_sink_constant_extent", m_extent);
        ddc::expose_to_pdi("krook_sink_constant_stiffness", m_stiffness);
        ddc::expose_to_pdi("krook_sink_constant_amplitude", m_amplitude);
        ddc::expose_to_pdi("krook_sink_constant_density", m_density);
        ddc::expose_to_pdi("krook_sink_constant_temperature", m_temperature);
        ddc::expose_to_pdi("krook_sink_constant_ftarget", ftarget_host);
        ddc::expose_to_pdi("krook_sink_constant_mask", mask_host);
        break;
    }
}

device_t<DSpanSpXVx> KrookSourceConstant::operator()(
        device_t<DSpanSpXVx> const allfdistribu,
        double const dt) const
{
    Kokkos::Profiling::pushRegion("KrookSource");

    auto ftarget = m_ftarget.span_view();
    auto mask = m_mask.span_view();
    auto const& amplitude = m_amplitude;

    ddc::for_each(
            ddc::policies::parallel_device,
            allfdistribu.domain(),
            KOKKOS_LAMBDA(IndexSpXVx const ispxvx) {
                allfdistribu(ispxvx)
                        = ftarget(ddc::select<IDimVx>(ispxvx))
                          + (allfdistribu(ispxvx) - ftarget(ddc::select<IDimVx>(ispxvx)))
                                    * Kokkos::exp(
                                            -amplitude * mask(ddc::select<IDimX>(ispxvx)) * dt);
            });

    Kokkos::Profiling::popRegion();
    return allfdistribu;
}
