#include <iomanip>

#include <fluid_moments.hpp>
#include <maxwellianequilibrium.hpp>
#include <pdi.h>
#include <rk2.hpp>

#include "collisions_inter.hpp"
#include "collisions_utils.hpp"

CollisionsInter::CollisionsInter(IDomainSpXVx const& mesh, double nustar0)
    : m_nustar0(nustar0)
    , m_nustar_profile_alloc(ddc::select<IDimSp, IDimX>(mesh))
{
    // validity checks
    if (ddc::select<IDimSp>(mesh).size() != 2) {
        throw std::runtime_error("Inter species collisions requires two kinetic species.");
    }
    if (m_nustar0 == 0.) {
        throw std::invalid_argument("Collision operator should not be used with nustar0=0.");
    }

    m_nustar_profile = m_nustar_profile_alloc.span_view();
    compute_nustar_profile(m_nustar_profile, m_nustar0);
    ddc::expose_to_pdi("collinter_nustar0", m_nustar0);
}

void CollisionsInter::get_derivative(
        device_t<DSpanSpXVx> const df,
        device_t<DViewSpXVx> const allfdistribu) const
{
    IDomainSpX grid_sp_x = allfdistribu.domain<IDimSp, IDimX>();
    device_t<DFieldSpX> density_f(grid_sp_x);
    device_t<DFieldSpX> fluid_velocity_f(grid_sp_x);
    device_t<DFieldSpX> temperature_f(grid_sp_x);
    auto density = density_f.span_view();
    auto fluid_velocity = fluid_velocity_f.span_view();
    auto temperature = temperature_f.span_view();

    DFieldVx const quadrature_coeffs(
            trapezoid_quadrature_coefficients(ddc::get_domain<IDimVx>(allfdistribu)));
    auto quadrature_coeffs_alloc = ddc::create_mirror_view_and_copy(
            Kokkos::DefaultExecutionSpace(),
            quadrature_coeffs.span_view());
    auto quadrature_coeffs_device = quadrature_coeffs_alloc.span_view();

    //Moments computation
    ddc::fill(density, 0.);
    ddc::for_each(
            ddc::policies::parallel_device,
            grid_sp_x,
            KOKKOS_LAMBDA(IndexSpX const ispx) {
                IndexSp isp(ddc::select<IDimSp>(ispx));
                IndexX ix(ddc::select<IDimX>(ispx));
                double particle_flux(0);
                double momentum_flux(0);
                for (IndexVx ivx : allfdistribu.domain<IDimVx>()) {
                    CoordVx const coordv = ddc::coordinate(ivx);
                    double const val(quadrature_coeffs_device(ivx) * allfdistribu(isp, ix, ivx));
                    density(isp, ix) += val;
                    particle_flux += val * coordv;
                    momentum_flux += val * coordv * coordv;
                }
                fluid_velocity(isp, ix) = particle_flux / density(isp, ix);
                temperature(isp, ix) = (momentum_flux - particle_flux * fluid_velocity(isp, ix))
                                       / density(isp, ix);
            });


    //Collision frequencies, momentum and energy exchange terms
    device_t<DFieldSpX> nustar_profile(grid_sp_x);
    ddc::deepcopy(nustar_profile, m_nustar_profile);
    device_t<DFieldSpX> collfreq_ab(grid_sp_x);
    device_t<DFieldSpX> momentum_exchange_ab_f(grid_sp_x);
    device_t<DFieldSpX> energy_exchange_ab_f(grid_sp_x);
    auto momentum_exchange_ab = momentum_exchange_ab_f.span_view();
    auto energy_exchange_ab = energy_exchange_ab_f.span_view();
    compute_collfreq_ab(collfreq_ab.span_view(), nustar_profile, density, temperature);
    compute_momentum_energy_exchange(
            momentum_exchange_ab.span_view(),
            energy_exchange_ab.span_view(),
            collfreq_ab,
            density,
            fluid_velocity,
            temperature);

    device_t<DFieldSpXVx> fmaxwellian_f(allfdistribu.domain());
    auto fmaxwellian = fmaxwellian_f.span_view();
    ddc::for_each(
            ddc::policies::parallel_device,
            allfdistribu.domain(),
            KOKKOS_LAMBDA(IndexSpXVx const ispxvx) {
                IndexSp isp(ddc::select<IDimSp>(ispxvx));
                IndexX ix(ddc::select<IDimX>(ispxvx));
                IndexVx ivx(ddc::select<IDimVx>(ispxvx));
                double const inv_sqrt_2piT = 1. / Kokkos::sqrt(2. * M_PI * temperature(isp, ix));
                CoordVx const vx = ddc::coordinate(ivx);
                fmaxwellian(ispxvx)
                        = density(isp, ix) * inv_sqrt_2piT
                          * Kokkos::exp(
                                  -(vx - fluid_velocity(isp, ix)) * (vx - fluid_velocity(isp, ix))
                                  / (2. * temperature(isp, ix)));
            });


    ddc::for_each(
            ddc::policies::parallel_device,
            allfdistribu.domain(),
            KOKKOS_LAMBDA(IndexSpXVx const ispxvx) {
                IndexSp isp(ddc::select<IDimSp>(ispxvx));
                IndexX ix(ddc::select<IDimX>(ispxvx));
                IndexVx ivx(ddc::select<IDimVx>(ispxvx));
                double const coordv = ddc::coordinate(ivx);
                double const term_v(coordv - fluid_velocity(isp, ix));
                df(isp, ix, ivx) = (2. * energy_exchange_ab(isp, ix)
                                            * (0.5 / temperature(isp, ix) * term_v * term_v - 0.5)
                                    + momentum_exchange_ab(isp, ix) * term_v)
                                   * fmaxwellian(isp, ix, ivx)
                                   / (density(isp, ix) * temperature(isp, ix));
            });
}


device_t<DSpanSpXVx> CollisionsInter::operator()(
        device_t<DSpanSpXVx> allfdistribu_device,
        double dt) const
{
    Kokkos::Profiling::pushRegion("CollisionsInter");
    RK2<device_t<DFieldSpXVx>> timestepper(allfdistribu_device.domain());

    timestepper
            .update(allfdistribu_device, dt, [&](device_t<DSpanSpXVx> dy, device_t<DViewSpXVx> y) {
                get_derivative(dy, y);
            });

    Kokkos::Profiling::popRegion();
    return allfdistribu_device;
}
