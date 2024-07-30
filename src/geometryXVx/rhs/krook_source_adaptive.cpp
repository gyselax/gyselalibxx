#include <cassert>
#include <stdexcept>
#include <string>

#include <ddc/ddc.hpp>

#include <maxwellianequilibrium.hpp>
#include <quadrature.hpp>
#include <rk2.hpp>
#include <species_info.hpp>
#include <trapezoid_quadrature.hpp>

#include "krook_source_adaptive.hpp"
#include "mask_tanh.hpp"


KrookSourceAdaptive::KrookSourceAdaptive(
        IDomainX const& gridx,
        IDomainVx const& gridvx,
        RhsType const type,
        double const extent,
        double const stiffness,
        double const amplitude,
        double const density,
        double const temperature)
    : m_type(type)
    , m_extent(extent)
    , m_stiffness(stiffness)
    , m_amplitude(amplitude)
    , m_density(density)
    , m_temperature(temperature)
    , m_mask(gridx)
    , m_ftarget(gridvx)
{
    // mask that defines the region where the operator is active
    host_t<DFieldX> mask_host(gridx);
    switch (m_type) {
    case RhsType::Source:
        // the mask equals one in the interval [x_left, x_right]
        mask_host = mask_tanh(gridx, m_extent, m_stiffness, MaskType::Normal, false);
        break;
    case RhsType::Sink:
        // the mask equals zero in the center of the plasma
        mask_host = mask_tanh(gridx, m_extent, m_stiffness, MaskType::Inverted, false);
        break;
    }
    ddc::parallel_deepcopy(m_mask.span_view(), mask_host);

    // target distribution function
    MaxwellianEquilibrium::compute_maxwellian(m_ftarget.span_view(), m_density, m_temperature, 0.);
    auto ftarget_host = ddc::create_mirror_view_and_copy(m_ftarget.span_view());

    switch (m_type) {
    case RhsType::Source:
        ddc::expose_to_pdi("krook_source_adaptive_extent", m_extent);
        ddc::expose_to_pdi("krook_source_adaptive_stiffness", m_stiffness);
        ddc::expose_to_pdi("krook_source_adaptive_amplitude", m_amplitude);
        ddc::expose_to_pdi("krook_source_adaptive_density", m_density);
        ddc::expose_to_pdi("krook_source_adaptive_temperature", m_temperature);
        ddc::expose_to_pdi("krook_source_adaptive_ftarget", ftarget_host);
        ddc::expose_to_pdi("krook_source_adaptive_mask", mask_host);
        break;
    case RhsType::Sink:
        ddc::expose_to_pdi("krook_sink_adaptive_extent", m_extent);
        ddc::expose_to_pdi("krook_sink_adaptive_stiffness", m_stiffness);
        ddc::expose_to_pdi("krook_sink_adaptive_amplitude", m_amplitude);
        ddc::expose_to_pdi("krook_sink_adaptive_density", m_density);
        ddc::expose_to_pdi("krook_sink_adaptive_temperature", m_temperature);
        ddc::expose_to_pdi("krook_sink_adaptive_ftarget", ftarget_host);
        ddc::expose_to_pdi("krook_sink_adaptive_mask", mask_host);
        break;
    }
}

void KrookSourceAdaptive::get_amplitudes(DSpanSpX amplitudes, DViewSpXVx const allfdistribu) const
{
    IdxRangeSp const dom_sp(ddc::get_domain<Species>(allfdistribu));
    assert(dom_sp.size() == 2);
    assert(charge(dom_sp.front()) * charge(dom_sp.back()) < 0.);
    std::optional<IdxSp> iion_opt;
    for (IdxSp const isp : dom_sp) {
        if (charge(isp) > 0.) {
            iion_opt = isp;
        }
    }
    IdxSp iion(iion_opt.value());
    IDomainVx const gridvx = allfdistribu.domain<IDimVx>();
    DFieldVx const quadrature_coeffs_alloc(
            trapezoid_quadrature_coefficients<Kokkos::DefaultExecutionSpace>(gridvx));
    auto quadrature_coeffs = quadrature_coeffs_alloc.span_view();

    auto const& amplitude = m_amplitude;
    auto const& density = m_density;
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            ddc::get_domain<IDimX>(allfdistribu),
            KOKKOS_LAMBDA(IndexX const ix) {
                amplitudes(IndexSpX(iion, ix)) = amplitude;
                double density_ion = 0.;
                double density_electron = 0.;

                for (IndexVx ivx : allfdistribu.domain<IDimVx>()) {
                    density_ion += quadrature_coeffs(ivx) * allfdistribu(iion, ix, ivx);
                    density_electron += quadrature_coeffs(ivx) * allfdistribu(ielec(), ix, ivx);
                }
                amplitudes(IndexSpX(ielec(), ix))
                        = amplitude * (density_ion - density) / (density_electron - density);
            });
}

void KrookSourceAdaptive::get_derivative(
        DSpanSpXVx df,
        DViewSpXVx allfdistribu,
        DViewSpXVx allfdistribu_start) const
{
    IDomainSpX grid_sp_x(allfdistribu.domain<Species, IDimX>());

    DFieldSpX amplitudes_alloc(grid_sp_x);
    auto amplitudes = amplitudes_alloc.span_view();
    get_amplitudes(amplitudes, allfdistribu);

    auto ftarget(m_ftarget.span_view());
    auto mask(m_mask.span_view());

    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            allfdistribu.domain(),
            KOKKOS_LAMBDA(IndexSpXVx const ispxvx) {
                IdxSp isp(ddc::select<Species>(ispxvx));
                IndexX ix(ddc::select<IDimX>(ispxvx));
                IndexVx ivx(ddc::select<IDimVx>(ispxvx));
                df(ispxvx) = -mask(ix) * amplitudes(isp, ix)
                             * (allfdistribu_start(ispxvx) - ftarget(ivx));
            });
}

DSpanSpXVx KrookSourceAdaptive::operator()(DSpanSpXVx const allfdistribu, double const dt) const
{
    Kokkos::Profiling::pushRegion("KrookSource");
    RK2<DFieldSpXVx> timestepper(allfdistribu.domain());

    timestepper.update(allfdistribu, dt, [&](DSpanSpXVx df, DViewSpXVx f) {
        get_derivative(df, f, allfdistribu);
    });
    Kokkos::Profiling::popRegion();
    return allfdistribu;
}
