#include <cassert>
#include <cmath>

#include <ddc/ddc.hpp>

#include <kinetic_source.hpp>
#include <species_info.hpp>

//TODO: this should be directly handled by ddc::Discretization really,
//      in the meantime, we do it ourselves
template <class IDim>
constexpr std::enable_if_t<!IDim::continuous_dimension_type::PERIODIC, double>
total_interval_length(DiscreteDomain<IDim> const& dom)
{
    return std::fabs(rlength(dom));
}

//TODO: this should be directly handled by ddc::Discretization really,
//      in the meantime, we do it ourselves
template <class RDim>
constexpr std::enable_if_t<RDim::PERIODIC, double> total_interval_length(
        DiscreteDomain<UniformPointSampling<RDim>> const& dom)
{
    return std::fabs(rlength(dom) + step<UniformPointSampling<RDim>>());
}

Kinetic_source::Kinetic_source(
        IDomainX const& gridx,
        IDomainVx const& gridvx,
        double const px_source,
        double const dx_source,
        double const source_amplitude,
        double const density_amplitude,
        double const energy_amplitude,
        double const temperature_source)
    : m_px_source(px_source)
    , m_dx_source(dx_source)
    , m_source_amplitude(source_amplitude)
    , m_density_amplitude(density_amplitude)
    , m_energy_amplitude(energy_amplitude)
    , m_temperature_source(temperature_source)
    , m_mask_source(gridx)
    , m_kinetic_source_vx(gridvx)
{
    // compute the source spatial extent (analoguous to a source mask)
    IVectX const Nx(gridx.size());
    CoordX const x_min(coordinate(gridx.front()));
    CoordX const x_max(coordinate(gridx.back()));
    double const Lx = total_interval_length(gridx);
    CoordX const x_scl(x_min + Lx * m_px_source);
    CoordX const x_scr(x_max - Lx * m_px_source);
    double const integral_coeff(
            m_dx_source
            * (std::log(std::cosh((Lx - x_scl) / m_dx_source))
               - std::log(std::cosh((Lx - x_scr) / m_dx_source))
               + std::log(std::cosh(x_scr / m_dx_source))
               - std::log(std::cosh(x_scl / m_dx_source))));
    for_each(policies::parallel_host, gridx, [=](IndexX const ix) {
        CoordX const coordx = coordinate(ix);
        m_mask_source(ix) = m_source_amplitude / integral_coeff
                            * (std::tanh((coordx - x_scl) / m_dx_source)
                               - std::tanh((coordx - x_scr) / m_dx_source));
    });

    // compute the source velocity profile (maxwellian profile here)
    double const coeff(std::sqrt(2 * M_PI * m_temperature_source));
    for_each(policies::parallel_host, gridvx, [=](IndexVx const ivx) {
        CoordVx const coordvx = coordinate(ivx);
        double const coordvx_sq = coordvx * coordvx;
        double const src_v_d_ivx = 1 / coeff * (1.5 - coordvx_sq / (2 * m_temperature_source))
                                   * std::exp(-coordvx_sq / (2 * m_temperature_source));
        double const src_v_e_ivx = -0.5 / coeff * (1 - coordvx_sq / m_temperature_source)
                                   * std::exp(-coordvx_sq / (2 * m_temperature_source));
        m_kinetic_source_vx(ivx)
                = m_density_amplitude * src_v_d_ivx + m_energy_amplitude * src_v_e_ivx;
    });
}

DSpanSpXVx Kinetic_source::operator()(DSpanSpXVx const allfdistribu, double const dt) const
{
    for_each(policies::parallel_host, allfdistribu.domain(), [=](IndexSpXVx const ispxvx) {
        double const df(
                m_mask_source(select<IDimX>(ispxvx)) * m_kinetic_source_vx(select<IDimVx>(ispxvx))
                * dt);
        allfdistribu(ispxvx) += df;
    });

    return allfdistribu;
}
