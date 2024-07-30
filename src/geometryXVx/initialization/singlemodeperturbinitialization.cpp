// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include "ddc_helper.hpp"
#include "singlemodeperturbinitialization.hpp"

SingleModePerturbInitialization::SingleModePerturbInitialization(
        DViewSpVx fequilibrium,
        host_t<IFieldSp> init_perturb_mode,
        host_t<DFieldSp> init_perturb_amplitude)
    : m_fequilibrium(fequilibrium)
    , m_init_perturb_mode(std::move(init_perturb_mode))
    , m_init_perturb_amplitude(std::move(init_perturb_amplitude))
{
}


DSpanSpXVx SingleModePerturbInitialization::operator()(DSpanSpXVx const allfdistribu) const
{
    IDomainSp const gridsp = allfdistribu.domain<IDimSp>();
    IDomainX const gridx = allfdistribu.domain<IDimX>();
    IDomainXVx const gridxvx = allfdistribu.domain<IDimX, IDimVx>();

    // Initialization of the perturbation
    DFieldX perturbation_alloc(gridx);
    ddc::ChunkSpan fequilibrium_proxy = m_fequilibrium.span_view();
    ddc::ChunkSpan perturbation_proxy = perturbation_alloc.span_view();
    ddc::for_each(gridsp, [&](IndexSp const isp) {
        perturbation_initialization(
                perturbation_proxy,
                m_init_perturb_mode(isp),
                m_init_perturb_amplitude(isp));
        // Initialization of the distribution function --> fill values
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                gridxvx,
                KOKKOS_LAMBDA(IndexXVx const ixvx) {
                    IndexX const ix = ddc::select<IDimX>(ixvx);
                    IndexVx const ivx = ddc::select<IDimVx>(ixvx);
                    double fdistribu_val
                            = fequilibrium_proxy(isp, ivx) * (1. + perturbation_proxy(ix));
                    if (fdistribu_val < 1.e-60) {
                        fdistribu_val = 1.e-60;
                    }
                    allfdistribu(isp, ix, ivx) = fdistribu_val;
                });
    });
    return allfdistribu;
}


SingleModePerturbInitialization SingleModePerturbInitialization::init_from_input(
        DViewSpVx allfequilibrium,
        IDomainSp dom_kinsp,
        PC_tree_t const& yaml_input_file)
{
    host_t<IFieldSp> init_perturb_mode(dom_kinsp);
    host_t<DFieldSp> init_perturb_amplitude(dom_kinsp);

    for (IndexSp const isp : dom_kinsp) {
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
        DSpanX const perturbation,
        int const perturb_mode,
        double const perturb_amplitude) const
{
    IDomainX const gridx = perturbation.domain();
    double const Lx = ddcHelper::total_interval_length(gridx);

    double const kx = perturb_mode * 2. * M_PI / Lx;
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            gridx,
            KOKKOS_LAMBDA(IndexX const ix) {
                CoordX const x = ddc::coordinate(ix);
                perturbation(ix) = perturb_amplitude * Kokkos::cos(kx * x);
            });
}
