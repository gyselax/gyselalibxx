// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include "ddc_helper.hpp"
#include "singlemodeperturbinitialization.hpp"


SingleModePerturbInitialization::SingleModePerturbInitialization(
        DViewSpVxVy fequilibrium,
        host_t<IFieldSp> init_perturb_mode,
        host_t<DFieldMemSp> init_perturb_amplitude)
    : m_fequilibrium(fequilibrium)
    , m_init_perturb_mode(std::move(init_perturb_mode))
    , m_init_perturb_amplitude(std::move(init_perturb_amplitude))
{
}

DSpanSpXYVxVy SingleModePerturbInitialization::operator()(DSpanSpXYVxVy const allfdistribu) const
{
    IdxRangeSp const gridsp = allfdistribu.domain<Species>();
    IDomainXY const gridxy = allfdistribu.domain<IDimX, IDimY>();
    IDomainXYVxVy const gridxyvxvy = allfdistribu.domain<IDimX, IDimY, IDimVx, IDimVy>();

    // Initialization of the perturbation
    DFieldXY perturbation_alloc(gridxy);
    ddc::ChunkSpan fequilibrium_proxy = m_fequilibrium.span_view();
    ddc::ChunkSpan perturbation_proxy = perturbation_alloc.span_view();
    ddc::for_each(gridsp, [&](IdxSp const isp) {
        perturbation_initialization(
                perturbation_proxy,
                m_init_perturb_mode(isp),
                m_init_perturb_amplitude(isp));

        // Initialization of the distribution function --> fill values
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                gridxyvxvy,
                KOKKOS_LAMBDA(IndexXYVxVy const ixyvxvy) {
                    IndexX const ix = ddc::select<IDimX>(ixyvxvy);
                    IndexY const iy = ddc::select<IDimY>(ixyvxvy);
                    IndexVx const ivx = ddc::select<IDimVx>(ixyvxvy);
                    IndexVy const ivy = ddc::select<IDimVy>(ixyvxvy);
                    double fdistribu_val
                            = fequilibrium_proxy(isp, ivx, ivy) * (1. + perturbation_proxy(ix, iy));
                    if (fdistribu_val < 1.e-60) {
                        fdistribu_val = 1.e-60;
                    }
                    allfdistribu(isp, ix, iy, ivx, ivy) = fdistribu_val;
                });
    });
    return allfdistribu;
}


SingleModePerturbInitialization SingleModePerturbInitialization::init_from_input(
        DViewSpVxVy allfequilibrium,
        IdxRangeSp dom_kinsp,
        PC_tree_t const& yaml_input_file)
{
    host_t<IFieldSp> init_perturb_mode(dom_kinsp);
    host_t<DFieldMemSp> init_perturb_amplitude(dom_kinsp);

    for (IdxSp const isp : dom_kinsp) {
        PC_tree_t const conf_isp = PCpp_get(yaml_input_file, ".SpeciesInfo[%d]", isp.uid());

        init_perturb_amplitude(isp) = PCpp_double(conf_isp, ".perturb_amplitude");
        init_perturb_mode(isp) = static_cast<int>(PCpp_int(conf_isp, ".perturb_mode"));
    }

    return SingleModePerturbInitialization(
            allfequilibrium,
            std::move(init_perturb_mode),
            std::move(init_perturb_amplitude));
}


void SingleModePerturbInitialization::perturbation_initialization(
        DSpanXY const perturbation,
        int const perturb_mode,
        double const perturb_amplitude) const
{
    IDomainXY const gridxy = perturbation.domain<IDimX, IDimY>();
    double const kx = perturb_mode * 2. * M_PI
                      / ddcHelper::total_interval_length(ddc::select<IDimX>(gridxy));
    double const ky = perturb_mode * 2. * M_PI
                      / ddcHelper::total_interval_length(ddc::select<IDimY>(gridxy));

    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            gridxy,
            KOKKOS_LAMBDA(IndexXY const ixy) {
                IndexX const ix = ddc::select<IDimX>(ixy);
                CoordX const x = ddc::coordinate(ix);
                IndexY const iy = ddc::select<IDimY>(ixy);
                CoordY const y = ddc::coordinate(iy);
                perturbation(ix, iy)
                        = perturb_amplitude * (Kokkos::cos(kx * x) + Kokkos::cos(ky * y));
            });
}
