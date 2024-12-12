#include <iomanip>

#include <ddc/ddc.hpp>
#include <ddc/pdi.hpp>

#include "collisions_inter.hpp"
#include "collisions_utils.hpp"
#include "fluid_moments.hpp"
#include "maxwellianequilibrium.hpp"
#include "rk2.hpp"

CollisionsInter::CollisionsInter(IdxRangeSpXVx const& mesh, double nustar0)
    : m_nustar0(nustar0)
    , m_nustar_profile_alloc(ddc::select<Species, GridX>(mesh))
{
    // validity checks
    if (ddc::select<Species>(mesh).size() != 2) {
        throw std::runtime_error("Inter species collisions requires two kinetic species.");
    }
    if (m_nustar0 == 0.) {
        throw std::invalid_argument("Collision operator should not be used with nustar0=0.");
    }

    m_nustar_profile = get_field(m_nustar_profile_alloc);
    compute_nustar_profile(m_nustar_profile, m_nustar0);
    ddc::expose_to_pdi("collinter_nustar0", m_nustar0);
}

void CollisionsInter::get_derivative(DFieldSpXVx const df, DConstFieldSpXVx const allfdistribu)
        const
{
    IdxRangeSpX grid_sp_x = get_idx_range<Species, GridX>(allfdistribu);
    DFieldMemSpX density_f(grid_sp_x);
    DFieldMemSpX fluid_velocity_f(grid_sp_x);
    DFieldMemSpX temperature_f(grid_sp_x);
    DFieldSpX density = get_field(density_f);
    DFieldSpX fluid_velocity = get_field(fluid_velocity_f);
    DFieldSpX temperature = get_field(temperature_f);

    IdxRangeVx const idx_range_vx(get_idx_range<GridVx>(allfdistribu));

    DFieldMemVx quadrature_coeffs_alloc(
            trapezoid_quadrature_coefficients<Kokkos::DefaultExecutionSpace>(idx_range_vx));
    DFieldVx quadrature_coeffs = get_field(quadrature_coeffs_alloc);

    //Moments computation
    ddc::parallel_fill(density, 0.);
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            grid_sp_x,
            KOKKOS_LAMBDA(IdxSpX const ispx) {
                IdxSp isp(ddc::select<Species>(ispx));
                IdxX ix(ddc::select<GridX>(ispx));
                double particle_flux(0);
                double momentum_flux(0);
                for (IdxVx ivx : idx_range_vx) {
                    CoordVx const coordv = ddc::coordinate(ivx);
                    double const val(quadrature_coeffs(ivx) * allfdistribu(isp, ix, ivx));
                    density(isp, ix) += val;
                    particle_flux += val * coordv;
                    momentum_flux += val * coordv * coordv;
                }
                fluid_velocity(isp, ix) = particle_flux / density(isp, ix);
                temperature(isp, ix) = (momentum_flux - particle_flux * fluid_velocity(isp, ix))
                                       / density(isp, ix);
            });


    //Collision frequencies, momentum and energy exchange terms
    DFieldMemSpX collfreq_ab(grid_sp_x);
    DFieldMemSpX momentum_exchange_ab_alloc(grid_sp_x);
    DFieldMemSpX energy_exchange_ab_alloc(grid_sp_x);
    DFieldSpX momentum_exchange_ab = get_field(momentum_exchange_ab_alloc);
    DFieldSpX energy_exchange_ab = get_field(energy_exchange_ab_alloc);
    compute_collfreq_ab(
            get_field(collfreq_ab),
            get_const_field(m_nustar_profile),
            get_const_field(density),
            get_const_field(temperature));
    compute_momentum_energy_exchange(
            momentum_exchange_ab,
            energy_exchange_ab,
            get_const_field(collfreq_ab),
            get_const_field(density),
            get_const_field(fluid_velocity),
            get_const_field(temperature));

    DFieldMemSpXVx fmaxwellian_alloc(get_idx_range(allfdistribu));
    DFieldSpXVx fmaxwellian = get_field(fmaxwellian_alloc);
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            get_idx_range(allfdistribu),
            KOKKOS_LAMBDA(IdxSpXVx const ispxvx) {
                IdxSp isp(ddc::select<Species>(ispxvx));
                IdxX ix(ddc::select<GridX>(ispxvx));
                IdxVx ivx(ddc::select<GridVx>(ispxvx));
                double const inv_sqrt_2piT = 1. / Kokkos::sqrt(2. * M_PI * temperature(isp, ix));
                CoordVx const vx = ddc::coordinate(ivx);
                fmaxwellian(ispxvx)
                        = density(isp, ix) * inv_sqrt_2piT
                          * Kokkos::exp(
                                  -(vx - fluid_velocity(isp, ix)) * (vx - fluid_velocity(isp, ix))
                                  / (2. * temperature(isp, ix)));
            });


    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            get_idx_range(allfdistribu),
            KOKKOS_LAMBDA(IdxSpXVx const ispxvx) {
                IdxSp isp(ddc::select<Species>(ispxvx));
                IdxX ix(ddc::select<GridX>(ispxvx));
                IdxVx ivx(ddc::select<GridVx>(ispxvx));
                double const coordv = ddc::coordinate(ivx);
                double const term_v(coordv - fluid_velocity(isp, ix));
                df(isp, ix, ivx) = (2. * energy_exchange_ab(isp, ix)
                                            * (0.5 / temperature(isp, ix) * term_v * term_v - 0.5)
                                    + momentum_exchange_ab(isp, ix) * term_v)
                                   * fmaxwellian(isp, ix, ivx)
                                   / (density(isp, ix) * temperature(isp, ix));
            });
}


DFieldSpXVx CollisionsInter::operator()(DFieldSpXVx allfdistribu, double dt) const
{
    Kokkos::Profiling::pushRegion("CollisionsInter");
    RK2<DFieldMemSpXVx> timestepper(get_idx_range(allfdistribu));

    timestepper.update(allfdistribu, dt, [&](DFieldSpXVx dy, DConstFieldSpXVx y) {
        get_derivative(dy, y);
    });

    Kokkos::Profiling::popRegion();
    return allfdistribu;
}
