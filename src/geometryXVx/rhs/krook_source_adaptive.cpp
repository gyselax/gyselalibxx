#include <cassert>
#include <stdexcept>
#include <string>

#include <ddc/ddc.hpp>

#include "krook_source_adaptive.hpp"
#include "mask_tanh.hpp"
#include "maxwellianequilibrium.hpp"
#include "quadrature.hpp"
#include "rk2.hpp"
#include "species_info.hpp"
#include "trapezoid_quadrature.hpp"


KrookSourceAdaptive::KrookSourceAdaptive(
        IdxRangeX const& gridx,
        IdxRangeVx const& gridvx,
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
    host_t<DFieldMemX> mask_host(gridx);
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
    ddc::parallel_deepcopy(get_field(m_mask), mask_host);

    // target distribution function
    MaxwellianEquilibrium::compute_maxwellian(get_field(m_ftarget), m_density, m_temperature, 0.);
    auto ftarget_host = ddc::create_mirror_view_and_copy(get_field(m_ftarget));

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

void KrookSourceAdaptive::get_amplitudes(DFieldSpX amplitudes, DConstFieldSpXVx const allfdistribu)
        const
{
    IdxRangeSp const dom_sp(get_idx_range<Species>(allfdistribu));
    assert(dom_sp.size() == 2);
    assert(charge(dom_sp.front()) * charge(dom_sp.back()) < 0.);
    std::optional<IdxSp> iion_opt;
    for (IdxSp const isp : dom_sp) {
        if (charge(isp) > 0.) {
            iion_opt = isp;
        }
    }
    IdxSp iion(iion_opt.value());
    IdxRangeVx const gridvx = get_idx_range<GridVx>(allfdistribu);
    DFieldMemVx const quadrature_coeffs_alloc(
            trapezoid_quadrature_coefficients<Kokkos::DefaultExecutionSpace>(gridvx));
    DConstFieldVx quadrature_coeffs = get_const_field(quadrature_coeffs_alloc);

    auto const& amplitude = m_amplitude;
    auto const& density = m_density;
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            get_idx_range<GridX>(allfdistribu),
            KOKKOS_LAMBDA(IdxX const ix) {
                amplitudes(IdxSpX(iion, ix)) = amplitude;
                double density_ion = 0.;
                double density_electron = 0.;

                for (IdxVx ivx : get_idx_range<GridVx>(allfdistribu)) {
                    density_ion += quadrature_coeffs(ivx) * allfdistribu(iion, ix, ivx);
                    density_electron += quadrature_coeffs(ivx) * allfdistribu(ielec(), ix, ivx);
                }
                amplitudes(IdxSpX(ielec(), ix))
                        = amplitude * (density_ion - density) / (density_electron - density);
            });
}

void KrookSourceAdaptive::get_derivative(
        DFieldSpXVx df,
        DConstFieldSpXVx allfdistribu,
        DConstFieldSpXVx allfdistribu_start) const
{
    IdxRangeSpX grid_sp_x(get_idx_range<Species, GridX>(allfdistribu));

    DFieldMemSpX amplitudes_alloc(grid_sp_x);
    auto amplitudes = get_field(amplitudes_alloc);
    get_amplitudes(amplitudes, allfdistribu);

    auto ftarget(get_field(m_ftarget));
    auto mask(get_field(m_mask));

    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            get_idx_range(allfdistribu),
            KOKKOS_LAMBDA(IdxSpXVx const ispxvx) {
                IdxSp isp(ddc::select<Species>(ispxvx));
                IdxX ix(ddc::select<GridX>(ispxvx));
                IdxVx ivx(ddc::select<GridVx>(ispxvx));
                df(ispxvx) = -mask(ix) * amplitudes(isp, ix)
                             * (allfdistribu_start(ispxvx) - ftarget(ivx));
            });
}

DFieldSpXVx KrookSourceAdaptive::operator()(DFieldSpXVx const allfdistribu, double const dt) const
{
    Kokkos::Profiling::pushRegion("KrookSource");
    RK2<DFieldMemSpXVx> timestepper(get_idx_range(allfdistribu));

    timestepper.update(allfdistribu, dt, [&](DFieldSpXVx df, DConstFieldSpXVx f) {
        get_derivative(df, f, allfdistribu);
    });
    Kokkos::Profiling::popRegion();
    return allfdistribu;
}
