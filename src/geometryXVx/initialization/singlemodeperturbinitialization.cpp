// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include "ddc_helper.hpp"
#include "singlemodeperturbinitialization.hpp"

SingleModePerturbInitialization::SingleModePerturbInitialization(
        DViewSpVx fequilibrium,
        host_t<ViewSp<int>> const init_perturb_mode,
        host_t<DViewSp> const init_perturb_amplitude)
    : m_fequilibrium(fequilibrium)
    , m_init_perturb_mode(init_perturb_mode)
    , m_init_perturb_amplitude(init_perturb_amplitude)
{
}

DSpanSpXVx SingleModePerturbInitialization::operator()(DSpanSpXVx const allfdistribu) const
{
    IDomainSp const gridsp = allfdistribu.domain<IDimSp>();
    IDomainX const gridx = allfdistribu.domain<IDimX>();
    IDomainXVx const gridxvx = allfdistribu.domain<IDimX, IDimVx>();

    // Initialization of the perturbation
    DFieldX perturbation_alloc(gridx);
    ddc::ChunkSpan perturbation = perturbation_alloc.span_view();
    ddc::for_each(gridsp, [&](IndexSp const isp) {
        perturbation_initialization(
                perturbation,
                m_init_perturb_mode(isp),
                m_init_perturb_amplitude(isp));
        // Initialization of the distribution function --> fill values
        ddc::for_each(
                ddc::policies::parallel_device,
                gridxvx,
                KOKKOS_CLASS_LAMBDA(IndexXVx const ixvx) {
                    IndexX const ix = ddc::select<IDimX>(ixvx);
                    IndexVx const ivx = ddc::select<IDimVx>(ixvx);
                    double fdistribu_val = m_fequilibrium(isp, ivx) * (1. + perturbation(ix));
                    if (fdistribu_val < 1.e-60) {
                        fdistribu_val = 1.e-60;
                    }
                    allfdistribu(isp, ix, ivx) = fdistribu_val;
                });
    });
    return allfdistribu;
}

void SingleModePerturbInitialization::perturbation_initialization(
        DSpanX const perturbation,
        int const mode,
        double const perturb_amplitude) const
{
    IDomainX const gridx = perturbation.domain();
    double const Lx = ddcHelper::total_interval_length(gridx);

    double const kx = mode * 2. * M_PI / Lx;
    ddc::for_each(
            ddc::policies::parallel_device,
            gridx,
            KOKKOS_LAMBDA(IndexX const ix) {
                CoordX const x = ddc::coordinate(ix);
                perturbation(ix) = perturb_amplitude * Kokkos::cos(kx * x);
            });
}
