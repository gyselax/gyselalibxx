// SPDX-License-Identifier: MIT

#include "diffusiveneutralsolver.hpp"
#include "quadrature.hpp"
#include "rk2.hpp"
#include "trapezoid_quadrature.hpp"


DiffusiveNeutralSolver::DiffusiveNeutralSolver(
        IReactionRate const& charge_exchange,
        IReactionRate const& ionization,
        IReactionRate const& recombination,
        double const normalization_coeff,
        SplineXBuilder_1d const& spline_x_builder,
        SplineXEvaluator_1d const& spline_x_evaluator,
        DConstFieldVx const& quadrature_coeffs)
    : m_charge_exchange(charge_exchange)
    , m_ionization(ionization)
    , m_recombination(recombination)
    , m_normalization_coeff(normalization_coeff)
    , m_spline_x_builder(spline_x_builder)
    , m_spline_x_evaluator(spline_x_evaluator)
    , m_quadrature_coeffs(quadrature_coeffs)
{
}

IdxSp DiffusiveNeutralSolver::find_ion(IdxRangeSp const dom_kinsp) const
{
    bool ion_found = false;
    IdxSp iion;
    for (IdxSp const isp : dom_kinsp) {
        if (charge(isp) > 0.) {
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

void DiffusiveNeutralSolver::get_derivative(
        DFieldSpMomX dn,
        DConstFieldSpMomX neutrals,
        DConstFieldSpX density,
        DConstFieldSpX velocity,
        DConstFieldSpX temperature) const
{
    IdxRangeSpX dom_fluidspx(get_idx_range<Species, GridX>(neutrals));

    // building reaction rates
    DFieldMemSpX charge_exchange_rate_alloc(dom_fluidspx);
    DFieldMemSpX ionization_rate_alloc(dom_fluidspx);
    DFieldMemSpX recombination_rate_alloc(dom_fluidspx);

    DFieldSpX charge_exchange_rate = get_field(charge_exchange_rate_alloc);
    DFieldSpX ionization_rate = get_field(ionization_rate_alloc);
    DFieldSpX recombination_rate = get_field(recombination_rate_alloc);

    m_charge_exchange(charge_exchange_rate, density, temperature);
    m_ionization(ionization_rate, density, temperature);
    m_recombination(recombination_rate, density, temperature);

    // compute diffusive model equation terms
    DFieldMemSpX density_equilibrium_velocity_alloc(dom_fluidspx);
    DFieldMemSpX diffusion_temperature_alloc(dom_fluidspx);
    DFieldSpX density_equilibrium_velocity = get_field(density_equilibrium_velocity_alloc);
    DFieldSpX diffusion_temperature = get_field(diffusion_temperature_alloc);

    IdxSp const iion(find_ion(get_idx_range<Species>(density)));
    IdxMom const ineutral_density(0);

    double const normalization_coeff_alpha0(m_normalization_coeff);
    double const mass_ratio(mass(ielec()) / mass(iion));
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            dom_fluidspx,
            KOKKOS_LAMBDA(IdxSpX const ifspx) {
                IdxSp const isp(ddc::select<Species>(ifspx));
                IdxX const ix(ddc::select<GridX>(ifspx));

                double const denom = density(iion, ix) * charge_exchange_rate(ifspx)
                                     + density(ielec(), ix) * ionization_rate(ifspx);

                density_equilibrium_velocity(ifspx)
                        = (density(ielec(), ix) * density(iion, ix) * recombination_rate(ifspx)
                           + neutrals(ifspx, ineutral_density) * density(iion, ix)
                                     * charge_exchange_rate(ifspx))
                          * velocity(iion, ix) * Kokkos::sqrt(mass_ratio) / denom;

                diffusion_temperature(ifspx)
                        = normalization_coeff_alpha0 * temperature(iion, ix) / (mass(isp) * denom);

                // density source is not solved here, we only solve transport.
            });

    // compute coordinates at which spatial derivatives are evaluated
    FieldMemX<CoordX> coords_eval_alloc(get_idx_range<GridX>(neutrals));
    auto coords_eval = get_field(coords_eval_alloc);
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            get_idx_range<GridX>(neutrals),
            KOKKOS_LAMBDA(IdxX const ix) { coords_eval(ix) = ddc::coordinate(ix); });

    // create chunks to store spatial derivatives
    DFieldMemSpX gradx_density_equilibrium_velocity_alloc(dom_fluidspx);
    DFieldMemSpX gradx_diffusion_temperature_alloc(dom_fluidspx);
    DFieldMemSpX gradx_neutrals_density_alloc(dom_fluidspx);
    DFieldMemSpX laplx_neutrals_density_alloc(dom_fluidspx);

    DFieldSpX gradx_density_equilibrium_velocity
            = get_field(gradx_density_equilibrium_velocity_alloc);
    DFieldSpX gradx_diffusion_temperature = get_field(gradx_diffusion_temperature_alloc);
    DFieldSpX gradx_neutrals_density = get_field(gradx_neutrals_density_alloc);
    DFieldSpX laplx_neutrals_density = get_field(laplx_neutrals_density_alloc);

    ddc::for_each(get_idx_range<Species>(neutrals), [&](IdxSp const isp) {
        // compute spline coefficients
        DBSFieldMemX density_equilibrium_velocity_spline_x_coeff(
                get_spline_idx_range(m_spline_x_builder));

        DBSFieldMemX diffusion_temperature_spline_x_coeff(get_spline_idx_range(m_spline_x_builder));

        DBSFieldMemX neutrals_density_spline_x_coeff(get_spline_idx_range(m_spline_x_builder));

        m_spline_x_builder(
                get_field(density_equilibrium_velocity_spline_x_coeff),
                get_const_field(density_equilibrium_velocity[isp]));

        m_spline_x_builder(
                get_field(diffusion_temperature_spline_x_coeff),
                get_const_field(diffusion_temperature[isp]));

        m_spline_x_builder(
                get_field(neutrals_density_spline_x_coeff),
                get_const_field(neutrals[IdxSpMom(isp, ineutral_density)]));

        // compute gradients
        m_spline_x_evaluator
                .deriv(get_field(gradx_density_equilibrium_velocity[isp]),
                       get_const_field(coords_eval),
                       get_const_field(density_equilibrium_velocity_spline_x_coeff));

        m_spline_x_evaluator
                .deriv(gradx_diffusion_temperature[isp],
                       get_const_field(coords_eval),
                       get_const_field(diffusion_temperature_spline_x_coeff));

        m_spline_x_evaluator
                .deriv(gradx_neutrals_density[isp],
                       get_const_field(coords_eval),
                       get_const_field(neutrals_density_spline_x_coeff));

        // compute laplacian
        DBSFieldMemX gradx_neutrals_density_spline_x_coeff(
                get_spline_idx_range(m_spline_x_builder));

        m_spline_x_builder(
                get_field(gradx_neutrals_density_spline_x_coeff),
                get_const_field(gradx_neutrals_density[isp]));

        m_spline_x_evaluator
                .deriv(laplx_neutrals_density[isp],
                       get_const_field(coords_eval),
                       get_const_field(gradx_neutrals_density_spline_x_coeff));
    });

    // build rhs of diffusive model equation
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            dom_fluidspx,
            KOKKOS_LAMBDA(IdxSpX const ifspx) {
                dn(ifspx, ineutral_density)
                        = -gradx_density_equilibrium_velocity(ifspx)
                          + gradx_diffusion_temperature(ifspx) * gradx_neutrals_density(ifspx)
                          + diffusion_temperature(ifspx) * laplx_neutrals_density(ifspx);
            }); // density source is not solved here, we only solve transport.
}

DFieldSpMomX DiffusiveNeutralSolver::operator()(
        DFieldSpMomX const neutrals,
        DConstFieldSpXVx const allfdistribu,
        DConstFieldX const efield,
        double const dt) const
{
    Kokkos::Profiling::pushRegion("DiffusiveNeutralSolver");
    RK2<DFieldMemSpMomX> timestepper(get_idx_range(neutrals));

    // moments computation
    IdxRangeSpX dom_kspx(get_idx_range(allfdistribu));
    DFieldMemSpX density_alloc(dom_kspx);
    DFieldMemSpX velocity_alloc(dom_kspx);
    DFieldMemSpX temperature_alloc(dom_kspx);

    DFieldSpX density = get_field(density_alloc);
    DFieldSpX velocity = get_field(velocity_alloc);
    DFieldSpX temperature = get_field(temperature_alloc);

    DConstFieldVx quadrature_coeffs = m_quadrature_coeffs;

    // fluid moments computation
    ddc::parallel_fill(density, 0.);
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            dom_kspx,
            KOKKOS_LAMBDA(IdxSpX const ispx) {
                double particle_flux(0);
                double momentum_flux(0);
                for (IdxVx const ivx : get_idx_range<GridVx>(allfdistribu)) {
                    CoordVx const coordv = ddc::coordinate(ivx);
                    double const val(quadrature_coeffs(ivx) * allfdistribu(ispx, ivx));
                    density(ispx) += val;
                    particle_flux += val * coordv;
                    momentum_flux += val * coordv * coordv;
                }
                velocity(ispx) = particle_flux / density(ispx);
                temperature(ispx)
                        = (momentum_flux - particle_flux * velocity(ispx)) / density(ispx);
            });

    timestepper.update(neutrals, dt, [&](DFieldSpMomX dn, DConstFieldSpMomX n) {
        get_derivative(dn, n, density, velocity, temperature);
    });
    Kokkos::Profiling::popRegion();
    return neutrals;
}
