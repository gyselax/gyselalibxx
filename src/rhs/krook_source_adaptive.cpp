#include <cassert>
#include <stdexcept>
#include <string>

#include <ddc/ddc.hpp>

#include <krook_source_adaptive.hpp>
#include <mask_tanh.hpp>
#include <maxwellianequilibrium.hpp>
#include <quadrature.hpp>
#include <trapezoid_quadrature.hpp>

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
                       double const time,
                       IndexSpX const ispx) { rhs(rhs_val, allfdistribu, time, ispx); });
        break;
    }

    // target distribution function
    MaxwellianEquilibrium::compute_maxwellian(m_ftarget, m_density, m_temperature, 0.);
}

double const KrookSourceAdaptive::get_amplitudes(DViewSpXVx allfdistribu, IndexSpX const ispx) const
{
    if (charge(select<IDimSp>(ispx)) >= 0.) {
        return m_amplitude;
    }

    assert(get_domain<IDimSp>(allfdistribu).size() == 2);
    IndexSp isp_ion;
    if (charge(select<IDimSp>(ispx)) != charge(get_domain<IDimSp>(allfdistribu).front())) {
        isp_ion = get_domain<IDimSp>(allfdistribu).front();
    } else {
        isp_ion = get_domain<IDimSp>(allfdistribu).back();
    }

    // compute densities for each species
    IDomainVx const gridvx = allfdistribu.domain<IDimVx>();
    Quadrature<IDimVx> const integrate_v(trapezoid_quadrature_coefficients(gridvx));
    double const density_ion = integrate_v(allfdistribu[isp_ion][select<IDimX>(ispx)]);
    double const density_electron = integrate_v(allfdistribu[select<IDimSp, IDimX>(ispx)]);

    double const amplitude
            = m_amplitude * (density_ion - m_density) / (density_electron - m_density);
    return amplitude;
}

void KrookSourceAdaptive::rhs(
        DSpanVx rhs,
        DViewSpXVx allfdistribu,
        double const time,
        IndexSpX const ispx) const
{
    double const amplitude = get_amplitudes(allfdistribu, ispx);
    for_each(policies::parallel_host, get_domain<IDimVx>(allfdistribu), [&](IndexVx const ivx) {
        IndexSpXVx ispxvx(select<IDimSp>(ispx), select<IDimX>(ispx), ivx);
        rhs(ivx) = -m_mask(select<IDimX>(ispx)) * amplitude
                   * (allfdistribu(ispxvx) - m_ftarget(ivx));
    });
}

DSpanSpXVx KrookSourceAdaptive::operator()(DSpanSpXVx allfdistribu, double const dt) const
{
    return (*m_solver)(allfdistribu, dt);
}
