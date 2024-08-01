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
        DConstFieldVx const& quadrature_coeffs) // for kinetic species
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

IdxSp KineticFluidCouplingSource::find_ion(IdxRangeSp const dom_kinsp) const
{
    bool ion_found = false;
    IdxSp iion;
    for (IdxSp const isp : dom_kinsp) {
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
        DFieldX density_source_neutral,
        DConstFieldSpX kinsp_density,
        DConstFieldSpMomX neutrals,
        DConstFieldSpX ionization,
        DConstFieldSpX recombination) const
{
    IdxSp const iion(find_ion(get_idx_range<Species>(kinsp_density)));
    // Neutrals density source computation
    IdxRangeSpMom const dom_msp(get_idx_range<Species, GridMom>(neutrals));
    IdxSpMom ineutral(dom_msp.front());
    IdxRangeSp const dom_sp(get_idx_range<Species>(neutrals));
    IdxSp ispneutral(dom_sp.front());

    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            get_idx_range(density_source_neutral),
            KOKKOS_LAMBDA(IdxX const ix) {
                density_source_neutral(ix) = neutrals(ineutral, ix) * kinsp_density(ielec(), ix)
                                                     * ionization(ispneutral, ix)
                                             - kinsp_density(iion, ix) * kinsp_density(ielec(), ix)
                                                       * recombination(ispneutral, ix);
            });
}

void KineticFluidCouplingSource::get_derivative_neutrals(
        DFieldSpMomX dn,
        DConstFieldSpMomX neutrals,
        DConstFieldX density_source_neutral) const
{
    // neutrals dn computation
    IdxRangeSpX dom_fluidspx(get_idx_range<Species, GridX>(neutrals));

    // compute diffusive model equation terms
    IdxMom const ineutral_density(0);
    double const normalization_coeff_alpha0(m_normalization_coeff);

    // build rhs of diffusive model equation
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            dom_fluidspx,
            KOKKOS_LAMBDA(IdxSpX const ifspx) {
                IdxX const ix(ifspx);
                dn(ifspx, ineutral_density)
                        = -density_source_neutral(ix) / normalization_coeff_alpha0;
            });
}

void KineticFluidCouplingSource::get_derivative_allfdistribu(
        DFieldSpXVx df,
        DConstFieldSpXVx allfdistribu,
        DConstFieldSpXVx velocity_shape_source) const
{
    // df computation
    IdxRangeSpXVx dom_kspx(get_idx_range(allfdistribu));

    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            dom_kspx,
            KOKKOS_LAMBDA(IdxSpXVx const ispxvx) { df(ispxvx) = velocity_shape_source(ispxvx); });
}

void KineticFluidCouplingSource::operator()(
        DFieldSpXVx const allfdistribu,
        DFieldSpMomX neutrals,
        double const dt) const
{
    Kokkos::Profiling::pushRegion("KineticFluidCouplingSource");
    RK2<DFieldMemSpMomX> timestepper_neutrals(get_idx_range(neutrals));
    RK2<DFieldMemSpXVx> timestepper_kinetic(get_idx_range(allfdistribu));

    // useful params and index ranges
    IdxSp const iion(find_ion(get_idx_range<Species>(allfdistribu)));

    // kinetic species fluid moments computation
    IdxRangeSpX dom_kspx(get_idx_range(allfdistribu));
    DFieldMemSpX kinsp_density_alloc(dom_kspx);
    DFieldMemSpX kinsp_velocity_alloc(dom_kspx);
    DFieldMemSpX kinsp_temperature_alloc(dom_kspx);

    DFieldSpX kinsp_density = get_field(kinsp_density_alloc);
    DFieldSpX kinsp_velocity = get_field(kinsp_velocity_alloc);
    DFieldSpX kinsp_temperature = get_field(kinsp_temperature_alloc);

    DConstFieldVx quadrature_coeffs = get_const_field(m_quadrature_coeffs);

    ddc::parallel_fill(Kokkos::DefaultExecutionSpace(), kinsp_density, 0.);
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            dom_kspx,
            KOKKOS_LAMBDA(IdxSpX const ispx) {
                double particle_flux(0);
                double momentum_flux(0);
                for (IdxVx const ivx : get_idx_range<GridVx>(allfdistribu)) {
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
    IdxRangeSpX dom_fluidspx(get_idx_range<Species, GridX>(neutrals));

    DFieldMemSpX ionization_rate_alloc(dom_fluidspx);
    DFieldMemSpX recombination_rate_alloc(dom_fluidspx);

    DFieldSpX ionization_rate = get_field(ionization_rate_alloc);
    DFieldSpX recombination_rate = get_field(recombination_rate_alloc);

    m_ionization(ionization_rate, kinsp_density, kinsp_temperature);
    m_recombination(recombination_rate, kinsp_density, kinsp_temperature);

    // source term computation
    IdxRangeX grid_x(get_idx_range<GridX>(allfdistribu));

    DFieldMemX density_source_neutral_alloc(grid_x);
    auto density_source_neutral = get_field(density_source_neutral_alloc);
    get_source_term(
            density_source_neutral,
            kinsp_density,
            neutrals,
            ionization_rate,
            recombination_rate);

    // S(v) velocity shape calculation for kinetic species
    DFieldMemSpXVx velocity_shape_source_alloc(get_idx_range(allfdistribu));
    DFieldSpXVx velocity_shape_source = get_field(velocity_shape_source_alloc);

    double density_coupling_coeff_proxy = m_density_coupling_coeff;
    double momentum_coupling_coeff_proxy = m_momentum_coupling_coeff;
    double energy_coupling_coeff_proxy = m_energy_coupling_coeff;
    double const normalization_coeff_alpha0_proxy = m_normalization_coeff;

    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            get_idx_range(allfdistribu),
            KOKKOS_LAMBDA(IdxSpXVx const ispxvx) {
                IdxX const ix(ispxvx);
                IdxVx const ivx(ispxvx);
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

    timestepper_kinetic.update(allfdistribu, dt, [&](DFieldSpXVx df, DConstFieldSpXVx f) {
        get_derivative_allfdistribu(df, f, velocity_shape_source);
    });
    timestepper_neutrals.update(neutrals, dt, [&](DFieldSpMomX dn, DConstFieldSpMomX n) {
        get_derivative_neutrals(dn, n, density_source_neutral);
    });

    Kokkos::Profiling::popRegion();
}
