#include <ddc/ddc.hpp>

#include <kinetic_source.hpp>
#include <mask_tanh.hpp>
#include <species_info.hpp>

/**
 * Solves the equation \partial f / \partial t = S(v, x),
 * where S is a space and velocity dependent source 
 * that does not depend on f or on time. 
 * Therefore f(t+dt) = f(t) + S*dt
 *
 * S = spatial_extent(x) * velocity_shape(v)
 * spatial_extent defines the location where the source is active.
 * spatial_extent is normalized, so that its integral along the spatial direction 
 * equals one. It has a hyperbolic tangent shape. It is equal to one in a central 
 * zone of the plasma of width defined by the extent parameter.
 *
 * velocity_shape defines the velocity profile of the source 
 * in the parallel velocity direction. It is the sum of a source that 
 * injects only density, and a source that injects only energy. If the density 
 * and energy parameters are equal to one (usual case), the resulting velocity_shape 
 * is maxwellian. 
 */
KineticSource::KineticSource(
        IDomainX const& gridx,
        IDomainVx const& gridvx,
        double const extent,
        double const stiffness,
        double const amplitude,
        double const density,
        double const energy,
        double const temperature)
    : m_extent(extent)
    , m_stiffness(stiffness)
    , m_amplitude(amplitude)
    , m_density(density)
    , m_energy(energy)
    , m_temperature(temperature)
    , m_spatial_extent(mask_tanh(gridx, extent, stiffness, MaskType::Normal, true))
    , m_velocity_shape(gridvx)
{
    // compute the source velocity profile (maxwellian profile if density = energy = 1.)
    double const coeff(1.0 / std::sqrt(2 * M_PI * m_temperature));
    for_each(policies::parallel_host, gridvx, [=](IndexVx const ivx) {
        CoordVx const coordvx = coordinate(ivx);
        double const coordvx_sq = coordvx * coordvx;
        double const density_source = coeff * (1.5 - coordvx_sq / (2 * m_temperature))
                                      * std::exp(-coordvx_sq / (2 * m_temperature));
        double const energy_source = -0.5 * coeff * (1 - coordvx_sq / m_temperature)
                                     * std::exp(-coordvx_sq / (2 * m_temperature));
        m_velocity_shape(ivx) = m_density * density_source + m_energy * energy_source;
    });
}

DSpanSpXVx KineticSource::operator()(DSpanSpXVx const allfdistribu, double const dt) const
{
    for_each(policies::parallel_host, allfdistribu.domain(), [=](IndexSpXVx const ispxvx) {
        double const df(
                m_amplitude * m_spatial_extent(select<IDimX>(ispxvx))
                * m_velocity_shape(select<IDimVx>(ispxvx)) * dt);
        allfdistribu(ispxvx) += df;
    });

    return allfdistribu;
}
