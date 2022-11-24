#include <ddc/ddc.hpp>

#include <kinetic_source.hpp>
#include <mask_tanh.hpp>
#include <species_info.hpp>

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
    , m_kinetic_source_vx(gridvx)
    , m_mask_source(mask_tanh(gridx, extent, stiffness, MaskType::Normal, true))
{
    // compute the source spatial extent (analoguous to a source mask)


    // compute the source velocity profile (maxwellian profile here)
    double const coeff(1.0 / std::sqrt(2 * M_PI * m_temperature));
    for_each(policies::parallel_host, gridvx, [=](IndexVx const ivx) {
        CoordVx const coordvx = coordinate(ivx);
        double const coordvx_sq = coordvx * coordvx;
        double const src_v_d_ivx = coeff * (1.5 - coordvx_sq / (2 * m_temperature))
                                   * std::exp(-coordvx_sq / (2 * m_temperature));
        double const src_v_e_ivx = -0.5 * coeff * (1 - coordvx_sq / m_temperature)
                                   * std::exp(-coordvx_sq / (2 * m_temperature));
        m_kinetic_source_vx(ivx) = m_density * src_v_d_ivx + m_energy * src_v_e_ivx;
    });
}

DSpanSpXVx KineticSource::operator()(DSpanSpXVx const allfdistribu, double const dt) const
{
    for_each(policies::parallel_host, allfdistribu.domain(), [=](IndexSpXVx const ispxvx) {
        double const df(
                m_amplitude * m_mask_source(select<IDimX>(ispxvx))
                * m_kinetic_source_vx(select<IDimVx>(ispxvx)) * dt);
        allfdistribu(ispxvx) += df;
    });

    return allfdistribu;
}
