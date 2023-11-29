#include <iomanip>

#include <fluid_moments.hpp>
#include <maxwellianequilibrium.hpp>
#include <pdi.h>
#include <rk2.hpp>

#include "collisions_inter.hpp"
#include "collisions_utils.hpp"

CollisionsInter::CollisionsInter(IDomainSpXVx const& mesh, double nustar0)
    : m_nustar0(nustar0)
    , m_nustar_profile(ddc::select<IDimSp, IDimX>(mesh))
{
    // validity checks
    if (ddc::select<IDimSp>(mesh).size() != 2) {
        throw std::runtime_error("Inter species collisions requires two kinetic species.");
    }
    if (m_nustar0 == 0.) {
        throw std::invalid_argument("Collision operator should not be used with nustar0=0.");
    }

    compute_nustar_profile(m_nustar_profile.span_view(), m_nustar0);
    ddc::expose_to_pdi("collinter_nustar0", m_nustar0);
}

void CollisionsInter::compute_rhs(DSpanSpXVx const rhs, DViewSpXVx const allfdistribu) const
{
    IDomainSp const dom_sp(ddc::get_domain<IDimSp>(allfdistribu));
    IDomainVx const gridvx(ddc::get_domain<IDimVx>(allfdistribu));

    DFieldVx const quadrature_coeffs = trapezoid_quadrature_coefficients(gridvx);
    Quadrature<IDimVx> integrate(quadrature_coeffs);
    FluidMoments moments(integrate);

    ddc::for_each(ddc::get_domain<IDimX>(allfdistribu), [&](IndexX const ix) {
        //Moments computation
        DFieldSp density(dom_sp);
        DFieldSp fluid_velocity(dom_sp);
        DFieldSp temperature(dom_sp);
        ddc::for_each(dom_sp, [&](IndexSp const isp) {
            IndexSpX ispx(isp, ix);
            moments(density(isp), allfdistribu[ispx], FluidMoments::s_density);
            moments(fluid_velocity(isp),
                    allfdistribu[ispx],
                    density(isp),
                    FluidMoments::s_velocity);
            moments(temperature(isp),
                    allfdistribu[ispx],
                    density(isp),
                    fluid_velocity(isp),
                    FluidMoments::s_temperature);
        });

        //Collision frequencies, momentum and energy exchange terms
        DFieldSp collfreq_ab(dom_sp);
        DFieldSp nustar_profile_copy(dom_sp);
        ddc::deepcopy(nustar_profile_copy, m_nustar_profile[ix]);
        compute_collfreq_ab(collfreq_ab.span_view(), nustar_profile_copy, density, temperature);
        DFieldSp momentum_exchange_ab(dom_sp);
        DFieldSp energy_exchange_ab(dom_sp);
        compute_momentum_energy_exchange(
                momentum_exchange_ab.span_view(),
                energy_exchange_ab.span_view(),
                collfreq_ab,
                density,
                fluid_velocity,
                temperature);

        ddc::for_each(dom_sp, [&](IndexSp const isp) {
            device_t<DFieldVx> fmaxwellian_device(gridvx);
            MaxwellianEquilibrium::compute_maxwellian(
                    fmaxwellian_device.span_view(),
                    density(isp),
                    temperature(isp),
                    fluid_velocity(isp));
            auto fmaxwellian = ddc::create_mirror_view_and_copy(fmaxwellian_device.span_view());


            ddc::for_each(gridvx, [&](IndexVx const ivx) {
                double const coordv = ddc::coordinate(ivx);
                double const term_v(coordv - fluid_velocity(isp));
                rhs(isp, ix, ivx) = (2. * energy_exchange_ab(isp)
                                             * (0.5 / temperature(isp) * term_v * term_v - 0.5)
                                     + momentum_exchange_ab(isp) * term_v)
                                    * fmaxwellian(ivx) / (density(isp) * temperature(isp));
            });
        });
    });
}


device_t<DSpanSpXVx> CollisionsInter::operator()(
        device_t<DSpanSpXVx> allfdistribu_device,
        double dt) const
{
    Kokkos::Profiling::pushRegion("CollisionsInter");
    auto allfdistribu_alloc = ddc::create_mirror_view_and_copy(allfdistribu_device);
    ddc::ChunkSpan allfdistribu = allfdistribu_alloc.span_view();
    RK2<DFieldSpXVx> timestepper(allfdistribu.domain());
    timestepper.update(allfdistribu, dt, [&](DSpanSpXVx dy, DViewSpXVx y) { compute_rhs(dy, y); });
    ddc::deepcopy(allfdistribu_device, allfdistribu);
    Kokkos::Profiling::popRegion();
    return allfdistribu_device;
}
