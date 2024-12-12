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
    , m_spatial_extent(gridx)
    , m_velocity_shape(gridvx)
{
    host_t<DFieldMemX> spatial_extent_host
            = mask_tanh(gridx, extent, stiffness, MaskType::Normal, true);
    ddc::parallel_deepcopy(get_field(m_spatial_extent), get_const_field(spatial_extent_host));
    // compute the source velocity profile (maxwellian profile if density = energy = 1.)
    double const coeff(1.0 / std::sqrt(2 * M_PI * m_temperature));
    host_t<DFieldMemVx> velocity_shape_host(gridvx);
    ddc::for_each(gridvx, [&](IdxVx const ivx) {
        CoordVx const coordvx = ddc::coordinate(ivx);
        double const coordvx_sq = coordvx * coordvx;
        double const density_source = coeff * (1.5 - coordvx_sq / (2 * temperature))
                                      * std::exp(-coordvx_sq / (2 * temperature));
        double const energy_source = -0.5 * coeff * (1 - coordvx_sq / temperature)
                                     * std::exp(-coordvx_sq / (2 * temperature));
        velocity_shape_host(ivx) = density * density_source + energy * energy_source;
    });
    ddc::parallel_deepcopy(get_field(m_velocity_shape), velocity_shape_host);
    ddc::expose_to_pdi("kinetic_source_extent", spatial_extent_host);
    ddc::expose_to_pdi("kinetic_source_stiffness", stiffness);
    ddc::expose_to_pdi("kinetic_source_amplitude", m_amplitude);
    ddc::expose_to_pdi("kinetic_source_density", m_density);
    ddc::expose_to_pdi("kinetic_source_energy", m_energy);
    ddc::expose_to_pdi("kinetic_source_temperature", m_temperature);
    ddc::expose_to_pdi("kinetic_source_velocity_shape", velocity_shape_host);
    ddc::expose_to_pdi("kinetic_source_spatial_extent", spatial_extent_host);
}

DFieldSpXVx KineticSource::operator()(DFieldSpXVx const allfdistribu, double const dt) const
{
    Kokkos::Profiling::pushRegion("KineticSource");
    DConstFieldVx velocity_shape = get_field(m_velocity_shape);

    DConstFieldX spatial_extent = get_field(m_spatial_extent);

    double const amplitude = m_amplitude;

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
