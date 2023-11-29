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

    // target distribution function
    device_t<DFieldVx> ftarget_device(gridvx);
    MaxwellianEquilibrium::compute_maxwellian(ftarget_device, m_density, m_temperature, 0.);
    auto ftarget = ddc::create_mirror_view_and_copy(ftarget_device.span_view());
    ddc::deepcopy(m_ftarget, ftarget);

    switch (m_type) {
    case RhsType::Source:
        ddc::expose_to_pdi("krook_source_adaptive_extent", m_extent);
        ddc::expose_to_pdi("krook_source_adaptive_stiffness", m_stiffness);
        ddc::expose_to_pdi("krook_source_adaptive_amplitude", m_amplitude);
        ddc::expose_to_pdi("krook_source_adaptive_density", m_density);
        ddc::expose_to_pdi("krook_source_adaptive_temperature", m_temperature);
        ddc::expose_to_pdi("krook_source_adaptive_ftarget", m_ftarget);
        ddc::expose_to_pdi("krook_source_adaptive_mask", m_mask);
        break;
    case RhsType::Sink:
        ddc::expose_to_pdi("krook_sink_adaptive_extent", m_extent);
        ddc::expose_to_pdi("krook_sink_adaptive_stiffness", m_stiffness);
        ddc::expose_to_pdi("krook_sink_adaptive_amplitude", m_amplitude);
        ddc::expose_to_pdi("krook_sink_adaptive_density", m_density);
        ddc::expose_to_pdi("krook_sink_adaptive_temperature", m_temperature);
        ddc::expose_to_pdi("krook_sink_adaptive_ftarget", m_ftarget);
        ddc::expose_to_pdi("krook_sink_adaptive_mask", m_mask);
        break;
    }
}

void KrookSourceAdaptive::get_amplitudes(
        device_t<DSpanSpX> amplitudes,
        device_t<DViewSpXVx> const allfdistribu) const
{
    IDomainSp const dom_sp(ddc::get_domain<IDimSp>(allfdistribu));
    assert(dom_sp.size() == 2);
    assert(charge(dom_sp.front()) * charge(dom_sp.back()) < 0);
    std::optional<IndexSp> iion_opt;
    for (IndexSp const isp : dom_sp) {
        if (charge(isp) > 0) {
            iion_opt = isp;
        }
    }
    IndexSp iion(iion_opt.value());
    IDomainVx const gridvx = allfdistribu.domain<IDimVx>();
    DFieldVx const quadrature_coeffs(trapezoid_quadrature_coefficients(gridvx));
    auto quadrature_coeffs_alloc = ddc::create_mirror_view_and_copy(
            Kokkos::DefaultExecutionSpace(),
            quadrature_coeffs.span_view());
    auto quadrature_coeffs_device = quadrature_coeffs_alloc.span_view();

    auto const& amplitude = m_amplitude;
    auto const& density = m_density;
    ddc::for_each(
            ddc::policies::parallel_device,
            ddc::get_domain<IDimX>(allfdistribu),
            DDC_LAMBDA(IndexX const ix) {
                amplitudes(IndexSpX(iion, ix)) = amplitude;
                double density_ion = 0.;
                double density_electron = 0.;

                for (IndexVx ivx : allfdistribu.domain<IDimVx>()) {
                    density_ion += quadrature_coeffs_device(ivx) * allfdistribu(iion, ix, ivx);
                    density_electron
                            += quadrature_coeffs_device(ivx) * allfdistribu(ielec(), ix, ivx);
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
    auto df_device_alloc
            = ddc::create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), df.span_view());
    auto df_device = df_device_alloc.span_view();

    auto allfdistribu_device_alloc = ddc::
            create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), allfdistribu.span_view());
    auto allfdistribu_device = allfdistribu_device_alloc.span_view();

    auto allfdistribu_start_device_alloc = ddc::create_mirror_view_and_copy(
            Kokkos::DefaultExecutionSpace(),
            allfdistribu_start.span_view());
    auto allfdistribu_start_device = allfdistribu_start_device_alloc.span_view();

    device_t<DFieldSpX> amplitudes_device_alloc(
            ddc::get_domain<IDimSp, IDimX>(allfdistribu_device));
    auto amplitudes_device = amplitudes_device_alloc.span_view();
    get_amplitudes(amplitudes_device, allfdistribu_device);

    auto mask_device_alloc
            = ddc::create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), m_mask.span_view());
    auto mask_device = mask_device_alloc.span_view();

    auto ftarget_device_alloc = ddc::
            create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), m_ftarget.span_view());
    auto ftarget_device = ftarget_device_alloc.span_view();

    ddc::for_each(
            ddc::policies::parallel_device,
            allfdistribu_device.domain(),
            DDC_LAMBDA(IndexSpXVx const ispxvx) {
                IndexSp isp(ddc::select<IDimSp>(ispxvx));
                IndexX ix(ddc::select<IDimX>(ispxvx));
                IndexVx ivx(ddc::select<IDimVx>(ispxvx));
                df_device(ispxvx) = -mask_device(ix) * amplitudes_device(isp, ix)
                                    * (allfdistribu_start_device(ispxvx) - ftarget_device(ivx));
            });
    ddc::deepcopy(df, df_device);
}

device_t<DSpanSpXVx> KrookSourceAdaptive::operator()(
        device_t<DSpanSpXVx> const allfdistribu_device,
        double const dt) const
{
    Kokkos::Profiling::pushRegion("KrookSource");
    auto allfdistribu_alloc = ddc::create_mirror_view_and_copy(allfdistribu_device);
    ddc::ChunkSpan allfdistribu = allfdistribu_alloc.span_view();
    RK2<DFieldSpXVx> timestepper(allfdistribu.domain());
    timestepper.update(allfdistribu, dt, [&](DSpanSpXVx df, DViewSpXVx f) {
        get_derivative(df, f, allfdistribu);
    });
    ddc::deepcopy(allfdistribu_device, allfdistribu);
    Kokkos::Profiling::popRegion();
    return allfdistribu_device;
}
