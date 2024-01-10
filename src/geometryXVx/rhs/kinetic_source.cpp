#include <ddc/ddc.hpp>

#include <species_info.hpp>

#include "kinetic_source.hpp"
#include "mask_tanh.hpp"

KineticSource::KineticSource(
        IDomainX const& gridx,
        IDomainVx const& gridvx,
        double const extent,
        double const stiffness,
        double const amplitude,
        double const density,
        double const energy,
        double const temperature)
    : m_amplitude(amplitude)
    , m_density(density)
    , m_energy(energy)
    , m_temperature(temperature)
    , m_spatial_extent(mask_tanh(gridx, extent, stiffness, MaskType::Normal, true))
    , m_velocity_shape(gridvx)
{
    // compute the source velocity profile (maxwellian profile if density = energy = 1.)
    double const coeff(1.0 / std::sqrt(2 * M_PI * m_temperature));
    ddc::for_each(gridvx, [=](IndexVx const ivx) {
        CoordVx const coordvx = ddc::coordinate(ivx);
        double const coordvx_sq = coordvx * coordvx;
        double const density_source = coeff * (1.5 - coordvx_sq / (2 * m_temperature))
                                      * std::exp(-coordvx_sq / (2 * m_temperature));
        double const energy_source = -0.5 * coeff * (1 - coordvx_sq / m_temperature)
                                     * std::exp(-coordvx_sq / (2 * m_temperature));
        m_velocity_shape(ivx) = m_density * density_source + m_energy * energy_source;
    });
    ddc::expose_to_pdi("kinetic_source_extent", m_spatial_extent);
    ddc::expose_to_pdi("kinetic_source_stiffness", stiffness);
    ddc::expose_to_pdi("kinetic_source_amplitude", m_amplitude);
    ddc::expose_to_pdi("kinetic_source_density", m_density);
    ddc::expose_to_pdi("kinetic_source_energy", m_energy);
    ddc::expose_to_pdi("kinetic_source_temperature", m_temperature);
    ddc::expose_to_pdi("kinetic_source_velocity_shape", m_velocity_shape);
    ddc::expose_to_pdi("kinetic_source_spatial_extent", m_spatial_extent);
}

device_t<DSpanSpXVx> KineticSource::operator()(
        device_t<DSpanSpXVx> const allfdistribu,
        double const dt) const
{
    Kokkos::Profiling::pushRegion("KineticSource");
    auto velocity_shape_device_alloc = ddc::create_mirror_view_and_copy(
            Kokkos::DefaultExecutionSpace(),
            m_velocity_shape.span_view());
    auto velocity_shape_device = velocity_shape_device_alloc.span_view();

    auto spatial_extent_device_alloc = ddc::create_mirror_view_and_copy(
            Kokkos::DefaultExecutionSpace(),
            m_spatial_extent.span_view());
    auto spatial_extent_device = spatial_extent_device_alloc.span_view();

    auto const& amplitude = m_amplitude;

    ddc::for_each(
            ddc::policies::parallel_device,
            allfdistribu.domain(),
            KOKKOS_LAMBDA(IndexSpXVx const ispxvx) {
                double const df(
                        amplitude * spatial_extent_device(ddc::select<IDimX>(ispxvx))
                        * velocity_shape_device(ddc::select<IDimVx>(ispxvx)) * dt);
                allfdistribu(ispxvx) += df;
            });

    Kokkos::Profiling::popRegion();
    return allfdistribu;
}
