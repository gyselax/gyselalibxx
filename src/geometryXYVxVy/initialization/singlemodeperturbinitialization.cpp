// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include "ddc_helper.hpp"
#include "singlemodeperturbinitialization.hpp"



SingleModePerturbInitialization::SingleModePerturbInitialization(
        DViewSpVxVy fequilibrium,
        ViewSp<int> const init_perturb_mode,
        DViewSp const init_perturb_amplitude)
    : m_fequilibrium(fequilibrium)
    , m_init_perturb_mode(init_perturb_mode)
    , m_init_perturb_amplitude(init_perturb_amplitude)
{
}

DSpanSpXYVxVy SingleModePerturbInitialization::operator()(DSpanSpXYVxVy const allfdistribu) const
{
    IDomainSp const gridsp = allfdistribu.domain<IDimSp>();
    IDomainXY const gridxy = allfdistribu.domain<IDimX, IDimY>();
    IDomainXYVxVy const gridxyvxvy = allfdistribu.domain<IDimX, IDimY, IDimVx, IDimVy>();

    // Initialization of the perturbation
    DFieldXY perturbation(gridxy);
    ddc::for_each(gridsp, [&](IndexSp const isp) {
        perturbation_initialization(
                perturbation,
                m_init_perturb_mode(isp),
                m_init_perturb_amplitude(isp));

        // Initialization of the distribution function --> fill values
        ddc::for_each(gridxyvxvy, [&](IndexXYVxVy const ixyvxvy) {
            IndexX const ix = ddc::select<IDimX>(ixyvxvy);
            IndexY const iy = ddc::select<IDimY>(ixyvxvy);
            IndexVx const ivx = ddc::select<IDimVx>(ixyvxvy);
            IndexVy const ivy = ddc::select<IDimVy>(ixyvxvy);
            double fdistribu_val = m_fequilibrium(isp, ivx, ivy) * (1. + perturbation(ix, iy));
            if (fdistribu_val < 1.e-60) {
                fdistribu_val = 1.e-60;
            }
            allfdistribu(isp, ix, iy, ivx, ivy) = fdistribu_val;
        });
    });
    return allfdistribu;
}

void SingleModePerturbInitialization::perturbation_initialization(
        DSpanXY const perturbation,
        int const mode,
        double const perturb_amplitude) const
{
    IDomainXY const gridxy = perturbation.domain<IDimX, IDimY>();
    double const kx
            = mode * 2. * M_PI / ddcHelper::total_interval_length(ddc::select<IDimX>(gridxy));
    double const ky
            = mode * 2. * M_PI / ddcHelper::total_interval_length(ddc::select<IDimY>(gridxy));

    ddc::for_each(gridxy, [&](IndexXY const ixy) {
        IndexX const ix = ddc::select<IDimX>(ixy);
        CoordX const x = ddc::coordinate(ix);
        IndexY const iy = ddc::select<IDimY>(ixy);
        CoordY const y = ddc::coordinate(iy);
        perturbation(ix, iy) = perturb_amplitude * (cos(kx * x) + cos(ky * y));
    });
}
