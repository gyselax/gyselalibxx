#include <ddc/ddc.hpp>

#include "kinetic_source.hpp"
#include "mask_tanh.hpp"
#include "species_info.hpp"

KineticSource::KineticSource(
        IdxRangeX const& gridx,
        IdxRangeVx const& gridvx,
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
    ddc::for_each(gridvx, [=](IdxVx const ivx) {
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

DFieldSpXVx KineticSource::operator()(DFieldSpXVx const allfdistribu, double const dt) const
{
    Kokkos::Profiling::pushRegion("KineticSource");
    auto velocity_shape_alloc = ddc::create_mirror_view_and_copy(
            Kokkos::DefaultExecutionSpace(),
            get_field(m_velocity_shape));
    auto velocity_shape = get_field(velocity_shape_alloc);

    auto spatial_extent_alloc = ddc::create_mirror_view_and_copy(
            Kokkos::DefaultExecutionSpace(),
            get_field(m_spatial_extent));
    auto spatial_extent = get_field(spatial_extent_alloc);

    auto const& amplitude = m_amplitude;

    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            get_idx_range(allfdistribu),
            KOKKOS_LAMBDA(IdxSpXVx const ispxvx) {
                double const df(
                        amplitude * spatial_extent(ddc::select<GridX>(ispxvx))
                        * velocity_shape(ddc::select<GridVx>(ispxvx)) * dt);
                allfdistribu(ispxvx) += df;
            });

    Kokkos::Profiling::popRegion();
    return allfdistribu;
}
