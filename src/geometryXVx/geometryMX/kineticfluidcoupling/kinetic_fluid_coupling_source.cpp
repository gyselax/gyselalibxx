#include <ddc/ddc.hpp>

#include <rk2.hpp>
#include <species_info.hpp>
#include <trapezoid_quadrature.hpp>

#include "kinetic_fluid_coupling_source.hpp"

KineticFluidCouplingSource::KineticFluidCouplingSource(
        double const density_coupling_coeff,
        double const momentum_coupling_coeff,
        double const energy_coupling_coeff,
        IReactionRate const& ionization,
        IReactionRate const& recombination,
        double const normalization_coeff,
        DViewVx const& quadrature_coeffs) // for kinetic species
    : m_density_coupling_coeff(density_coupling_coeff)
    , m_momentum_coupling_coeff(momentum_coupling_coeff)
    , m_energy_coupling_coeff(energy_coupling_coeff)
    , m_ionization(ionization)
    , m_recombination(recombination)
    , m_normalization_coeff(normalization_coeff)
    , m_quadrature_coeffs(quadrature_coeffs)
{
    ddc::expose_to_pdi(
            "kinetic_fluid_coupling_source_density_coupling_coeff",
            m_density_coupling_coeff);
    ddc::expose_to_pdi(
            "kinetic_fluid_coupling_source_momentum_coupling_coeff",
            m_momentum_coupling_coeff);
    ddc::expose_to_pdi(
            "kinetic_fluid_coupling_source_energy_coupling_coeff",
            m_energy_coupling_coeff);
}

IndexSp KineticFluidCouplingSource::find_ion(IDomainSp const dom_kinsp) const
{
    bool ion_found = false;
    IndexSp iion;
    for (IndexSp const isp : dom_kinsp) {
        if (charge(isp) > 0) {
            ion_found = true;
            iion = isp;
        }
    }
    if (!ion_found) {
        throw std::runtime_error("ion not found");
    }
    assert(dom_kinsp.size() == 2);

    return iion;
}

void KineticFluidCouplingSource::get_source_term(
        DSpanX density_source_neutral,
        DViewSpX kinsp_density,
        DViewSpMX neutrals,
        DViewSpX ionization,
        DViewSpX recombination) const
{
    IndexSp const iion(find_ion(ddc::get_domain<IDimSp>(kinsp_density)));
    // Neutrals density source computation
    IDomainSpM const dom_msp(ddc::get_domain<IDimSp, IDimM>(neutrals));
    IndexSpM ineutral(dom_msp.front());
    IDomainSp const dom_sp(ddc::get_domain<IDimSp>(neutrals));
    IndexSp ispneutral(dom_sp.front());

    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            density_source_neutral.domain(),
            KOKKOS_LAMBDA(IndexX const ix) {
                density_source_neutral(ix) = neutrals(ineutral, ix) * kinsp_density(ielec(), ix)
                                                     * ionization(ispneutral, ix)
                                             - kinsp_density(iion, ix) * kinsp_density(ielec(), ix)
                                                       * recombination(ispneutral, ix);
            });
}

void KineticFluidCouplingSource::get_derivative_neutrals(
        DSpanSpMX dn,
        DViewSpMX neutrals,
        DViewX density_source_neutral) const
{
    // neutrals dn computation
    IDomainSpX dom_fluidspx(ddc::get_domain<IDimSp, IDimX>(neutrals));

    // compute diffusive model equation terms
    IndexM const ineutral_density(0);
    double const normalization_coeff_alpha0(m_normalization_coeff);

    // build rhs of diffusive model equation
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            dom_fluidspx,
            KOKKOS_LAMBDA(IndexSpX const ifspx) {
                IndexX const ix(ifspx);
                dn(ifspx, ineutral_density)
                        = -density_source_neutral(ix) / normalization_coeff_alpha0;
            });
}

void KineticFluidCouplingSource::get_derivative_allfdistribu(
        DSpanSpXVx df,
        DViewSpXVx allfdistribu,
        DViewSpXVx velocity_shape_source) const
{
    // df computation
    IDomainSpXVx dom_kspx(allfdistribu.domain());

    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            dom_kspx,
            KOKKOS_LAMBDA(IndexSpXVx const ispxvx) { df(ispxvx) = velocity_shape_source(ispxvx); });
}

void KineticFluidCouplingSource::operator()(
        DSpanSpXVx const allfdistribu,
        DSpanSpMX neutrals,
        double const dt) const
{
    Kokkos::Profiling::pushRegion("KineticFluidCouplingSource");
    RK2<DFieldSpMX> timestepper_neutrals(neutrals.domain());
    RK2<DFieldSpXVx> timestepper_kinetic(allfdistribu.domain());

    // useful params and domains
    IndexSp const iion(find_ion(ddc::get_domain<IDimSp>(allfdistribu)));

    // kinetic species fluid moments computation
    IDomainSpX dom_kspx(allfdistribu.domain());
    DFieldSpX kinsp_density_alloc(dom_kspx);
    DFieldSpX kinsp_velocity_alloc(dom_kspx);
    DFieldSpX kinsp_temperature_alloc(dom_kspx);

    DSpanSpX kinsp_density = kinsp_density_alloc.span_view();
    DSpanSpX kinsp_velocity = kinsp_velocity_alloc.span_view();
    DSpanSpX kinsp_temperature = kinsp_temperature_alloc.span_view();

    DViewVx quadrature_coeffs = m_quadrature_coeffs.span_view();

    ddc::parallel_fill(Kokkos::DefaultExecutionSpace(), kinsp_density, 0.);
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            dom_kspx,
            KOKKOS_LAMBDA(IndexSpX const ispx) {
                double particle_flux(0);
                double momentum_flux(0);
                for (IndexVx const ivx : allfdistribu.domain<IDimVx>()) {
                    CoordVx const coordv = ddc::coordinate(ivx);
                    double const val(quadrature_coeffs(ivx) * allfdistribu(ispx, ivx));
                    kinsp_density(ispx) += val;
                    particle_flux += val * coordv;
                    momentum_flux += val * coordv * coordv;
                }
                kinsp_velocity(ispx) = particle_flux / kinsp_density(ispx);
                kinsp_temperature(ispx) = (momentum_flux - particle_flux * kinsp_velocity(ispx))
                                          / kinsp_density(ispx);
            });

    // building reaction rates
    IDomainSpX dom_fluidspx(ddc::get_domain<IDimSp, IDimX>(neutrals));

    DFieldSpX ionization_rate_alloc(dom_fluidspx);
    DFieldSpX recombination_rate_alloc(dom_fluidspx);

    DSpanSpX ionization_rate = ionization_rate_alloc.span_view();
    DSpanSpX recombination_rate = recombination_rate_alloc.span_view();

    m_ionization(ionization_rate, kinsp_density, kinsp_temperature);
    m_recombination(recombination_rate, kinsp_density, kinsp_temperature);

    // source term computation
    IDomainX grid_x(allfdistribu.domain<IDimX>());

    DFieldX density_source_neutral_alloc(grid_x);
    auto density_source_neutral = density_source_neutral_alloc.span_view();
    get_source_term(
            density_source_neutral,
            kinsp_density,
            neutrals,
            ionization_rate,
            recombination_rate);

    // S(v) velocity shape calculation for kinetic species
    DFieldSpXVx velocity_shape_source_alloc(allfdistribu.domain());
    DSpanSpXVx velocity_shape_source = velocity_shape_source_alloc.span_view();

    double density_coupling_coeff_proxy = m_density_coupling_coeff;
    double momentum_coupling_coeff_proxy = m_momentum_coupling_coeff;
    double energy_coupling_coeff_proxy = m_energy_coupling_coeff;
    double const normalization_coeff_alpha0_proxy = m_normalization_coeff;

    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            allfdistribu.domain(),
            KOKKOS_LAMBDA(IndexSpXVx const ispxvx) {
                IndexX const ix(ispxvx);
                IndexVx const ivx(ispxvx);
                CoordVx const coordvx = ddc::coordinate(ivx);
                double const neutral_temperature
                        = (kinsp_temperature(iion, ix) + kinsp_temperature(ielec(), ix)) / 2.;

                double const coordvx_sq = coordvx * coordvx;
                double const density_source
                        = density_coupling_coeff_proxy
                          * (1.5 - coordvx_sq / (2 * neutral_temperature))
                          * Kokkos::exp(-coordvx_sq / (2 * neutral_temperature));
                double const momentum_source
                        = momentum_coupling_coeff_proxy * Kokkos::sqrt(2 / neutral_temperature)
                          * coordvx * Kokkos::exp(-coordvx_sq / (2 * neutral_temperature));
                double const energy_source = 2 * energy_coupling_coeff_proxy
                                             * (-1 + coordvx_sq / neutral_temperature)
                                             * Kokkos::exp(-coordvx_sq / (2 * neutral_temperature));
                velocity_shape_source(ispxvx) = (density_source_neutral(ix)
                                                 / (Kokkos::sqrt(2 * M_PI * neutral_temperature)
                                                    * normalization_coeff_alpha0_proxy))
                                                        * density_source
                                                + momentum_source + energy_source;
            });

    timestepper_kinetic.update(allfdistribu, dt, [&](DSpanSpXVx df, DViewSpXVx f) {
        get_derivative_allfdistribu(df, f, velocity_shape_source);
    });
    timestepper_neutrals.update(neutrals, dt, [&](DSpanSpMX dn, DViewSpMX n) {
        get_derivative_neutrals(dn, n, density_source_neutral);
    });

    Kokkos::Profiling::popRegion();
}
