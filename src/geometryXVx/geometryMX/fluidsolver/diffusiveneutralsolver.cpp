// SPDX-License-Identifier: MIT

#include <quadrature.hpp>
#include <rk2.hpp>
#include <trapezoid_quadrature.hpp>

#include "diffusiveneutralsolver.hpp"


DiffusiveNeutralSolver::DiffusiveNeutralSolver(
        IReactionRate const& chargexchange,
        IReactionRate const& ionization,
        IReactionRate const& recombination,
        double const temperature,
        double const normalization_coeff,
        SplineXBuilder_1d const& spline_x_builder,
        SplineXEvaluator_1d const& spline_x_evaluator,
        DViewVx const& quadrature_coeffs)
    : m_chargexchange(chargexchange)
    , m_ionization(ionization)
    , m_recombination(recombination)
    , m_temperature(temperature)
    , m_normalization_coeff(normalization_coeff)
    , m_spline_x_builder(spline_x_builder)
    , m_spline_x_evaluator(spline_x_evaluator)
    , m_quadrature_coeffs(quadrature_coeffs)
{
}

IndexSp DiffusiveNeutralSolver::find_ion(IDomainSp const dom_kinsp) const
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

void DiffusiveNeutralSolver::get_derivative(
        DSpanSpMX dn,
        DViewSpMX neutrals,
        DViewSpX density,
        DViewSpX velocity,
        DViewSpX temperature) const
{
    IDomainSpX dom_fluidspx(ddc::get_domain<IDimSp, IDimX>(neutrals));

    // building reaction rates
    DFieldSpX chargexchange_rate_alloc(dom_fluidspx);
    DFieldSpX ionization_rate_alloc(dom_fluidspx);
    DFieldSpX recombination_rate_alloc(dom_fluidspx);

    DSpanSpX chargexchange_rate = chargexchange_rate_alloc.span_view();
    DSpanSpX ionization_rate = ionization_rate_alloc.span_view();
    DSpanSpX recombination_rate = recombination_rate_alloc.span_view();

    m_chargexchange(chargexchange_rate, density, temperature);
    m_ionization(ionization_rate, density, temperature);
    m_recombination(recombination_rate, density, temperature);

    // compute diffusive model equation terms
    DFieldSpX density_equilibrium_velocity_alloc(dom_fluidspx);
    DFieldSpX diffusion_temperature_alloc(dom_fluidspx);
    DFieldSpX density_source_alloc(dom_fluidspx);
    ddc::ChunkSpan density_equilibrium_velocity = density_equilibrium_velocity_alloc.span_view();
    ddc::ChunkSpan diffusion_temperature = diffusion_temperature_alloc.span_view();
    ddc::ChunkSpan density_source = density_source_alloc.span_view();

    IndexSp const iion(find_ion(ddc::get_domain<IDimSp>(density)));
    IndexM const ineutral_density(0);

    double const normalization_coeff_alpha0(m_normalization_coeff);
    double const temperature_neutrals(m_temperature);
    double const mass_ratio(mass(ielec()) / mass(iion));
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            dom_fluidspx,
            KOKKOS_LAMBDA(IndexSpX const ifspx) {
                IndexSp const isp(ddc::select<IDimSp>(ifspx));
                IndexX const ix(ddc::select<IDimX>(ifspx));

                double const denom = density(iion, ix) * chargexchange_rate(ifspx)
                                     + density(ielec(), ix) * ionization_rate(ifspx);

                density_equilibrium_velocity(ifspx)
                        = (density(ielec(), ix) * density(iion, ix) * recombination_rate(ifspx)
                           + neutrals(ifspx, ineutral_density) * density(iion, ix)
                                     * chargexchange_rate(ifspx))
                          * velocity(iion, ix) * Kokkos::sqrt(mass_ratio) / denom;

                diffusion_temperature(ifspx)
                        = normalization_coeff_alpha0 * temperature_neutrals / (mass(isp) * denom);

                density_source(ifspx)
                        = (density(ielec(), ix) * density(iion, ix) * recombination_rate(ifspx)
                           - neutrals(ifspx, ineutral_density) * density(ielec(), ix)
                                     * ionization_rate(ifspx))
                          / normalization_coeff_alpha0;
            });

    // compute coordinates at which spatial derivatives are evaluated
    FieldX<CoordX> coords_eval_alloc(neutrals.domain<IDimX>());
    auto coords_eval = coords_eval_alloc.span_view();
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            neutrals.domain<IDimX>(),
            KOKKOS_LAMBDA(IndexX const ix) { coords_eval(ix) = ddc::coordinate(ix); });

    // create chunks to store spatial derivatives
    DFieldSpX gradx_density_equilibrium_velocity_alloc(dom_fluidspx);
    DFieldSpX gradx_diffusion_temperature_alloc(dom_fluidspx);
    DFieldSpX gradx_neutrals_density_alloc(dom_fluidspx);
    DFieldSpX laplx_neutrals_density_alloc(dom_fluidspx);

    ddc::ChunkSpan gradx_density_equilibrium_velocity
            = gradx_density_equilibrium_velocity_alloc.span_view();
    ddc::ChunkSpan gradx_diffusion_temperature = gradx_diffusion_temperature_alloc.span_view();
    ddc::ChunkSpan gradx_neutrals_density = gradx_neutrals_density_alloc.span_view();
    ddc::ChunkSpan laplx_neutrals_density = laplx_neutrals_density_alloc.span_view();

    ddc::for_each(neutrals.domain<IDimSp>(), [&](IndexSp const isp) {
        // compute spline coefficients
        DBSFieldX density_equilibrium_velocity_spline_x_coeff(m_spline_x_builder.spline_domain());

        DBSFieldX diffusion_temperature_spline_x_coeff(m_spline_x_builder.spline_domain());

        DBSFieldX neutrals_density_spline_x_coeff(m_spline_x_builder.spline_domain());

        m_spline_x_builder(
                density_equilibrium_velocity_spline_x_coeff.span_view(),
                density_equilibrium_velocity[isp].span_cview());

        m_spline_x_builder(
                diffusion_temperature_spline_x_coeff.span_view(),
                diffusion_temperature[isp].span_cview());

        m_spline_x_builder(
                neutrals_density_spline_x_coeff.span_view(),
                neutrals[IndexSpM(isp, ineutral_density)].span_cview());

        // compute gradients
        m_spline_x_evaluator
                .deriv(gradx_density_equilibrium_velocity[isp].span_view(),
                       coords_eval.span_cview(),
                       density_equilibrium_velocity_spline_x_coeff.span_cview());

        m_spline_x_evaluator
                .deriv(gradx_diffusion_temperature[isp],
                       coords_eval.span_cview(),
                       diffusion_temperature_spline_x_coeff.span_cview());

        m_spline_x_evaluator
                .deriv(gradx_neutrals_density[isp],
                       coords_eval.span_cview(),
                       neutrals_density_spline_x_coeff.span_cview());

        // compute laplacian
        DBSFieldX gradx_neutrals_density_spline_x_coeff(m_spline_x_builder.spline_domain());

        m_spline_x_builder(
                gradx_neutrals_density_spline_x_coeff.span_view(),
                gradx_neutrals_density[isp].span_cview());

        m_spline_x_evaluator
                .deriv(laplx_neutrals_density[isp],
                       coords_eval.span_cview(),
                       gradx_neutrals_density_spline_x_coeff.span_cview());
    });

    // build rhs of diffusive model equation
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            dom_fluidspx,
            KOKKOS_LAMBDA(IndexSpX const ifspx) {
                dn(ifspx, ineutral_density)
                        = density_source(ifspx) - gradx_density_equilibrium_velocity(ifspx)
                          + gradx_diffusion_temperature(ifspx) * gradx_neutrals_density(ifspx)
                          + diffusion_temperature(ifspx) * laplx_neutrals_density(ifspx);
            });
}

DSpanSpMX DiffusiveNeutralSolver::operator()(
        DSpanSpMX const neutrals,
        DViewSpXVx const allfdistribu,
        DViewX const efield,
        double const dt) const
{
    RK2<DFieldSpMX> timestepper(neutrals.domain());

    // moments computation
    IDomainSpX dom_kspx(allfdistribu.domain());
    DFieldSpX density_alloc(dom_kspx);
    DFieldSpX velocity_alloc(dom_kspx);
    DFieldSpX temperature_alloc(dom_kspx);

    DSpanSpX density = density_alloc.span_view();
    DSpanSpX velocity = velocity_alloc.span_view();
    DSpanSpX temperature = temperature_alloc.span_view();

    DViewVx quadrature_coeffs = m_quadrature_coeffs.span_view();

    // fluid moments computation
    ddc::parallel_fill(density, 0.);
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            dom_kspx,
            KOKKOS_LAMBDA(IndexSpX const ispx) {
                double particle_flux(0);
                double momentum_flux(0);
                for (IndexVx const ivx : allfdistribu.domain<IDimVx>()) {
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

    timestepper.update(neutrals, dt, [&](DSpanSpMX dn, DViewSpMX n) {
        get_derivative(dn, n, density, velocity, temperature);
    });
    return neutrals;
}
