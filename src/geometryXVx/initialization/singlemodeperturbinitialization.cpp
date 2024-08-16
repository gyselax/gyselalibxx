// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include "ddc_helper.hpp"
#include "singlemodeperturbinitialization.hpp"

SingleModePerturbInitialization::SingleModePerturbInitialization(
        DConstFieldSpVx fequilibrium,
        host_t<IFieldMemSp> init_perturb_mode,
        host_t<DFieldMemSp> init_perturb_amplitude)
    : m_fequilibrium(fequilibrium)
    , m_init_perturb_mode(std::move(init_perturb_mode))
    , m_init_perturb_amplitude(std::move(init_perturb_amplitude))
{
}


DFieldSpXVx SingleModePerturbInitialization::operator()(DFieldSpXVx const allfdistribu) const
{
    IdxRangeSp const gridsp = get_idx_range<Species>(allfdistribu);
    IdxRangeX const gridx = get_idx_range<GridX>(allfdistribu);
    IdxRangeXVx const gridxvx = get_idx_range<GridX, GridVx>(allfdistribu);

    // Initialization of the perturbation
    DFieldMemX perturbation_alloc(gridx);
    DConstFieldSpVx fequilibrium_proxy = get_const_field(m_fequilibrium);
    DFieldX perturbation_proxy = get_field(perturbation_alloc);
    ddc::for_each(gridsp, [&](IdxSp const isp) {
        perturbation_initialization(
                perturbation_proxy,
                m_init_perturb_mode(isp),
                m_init_perturb_amplitude(isp));
        // Initialization of the distribution function --> fill values
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                gridxvx,
                KOKKOS_LAMBDA(IdxXVx const ixvx) {
                    IdxX const ix = ddc::select<GridX>(ixvx);
                    IdxVx const ivx = ddc::select<GridVx>(ixvx);
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
        DConstFieldSpVx allfequilibrium,
        IdxRangeSp idx_range_kinsp,
        PC_tree_t const& yaml_input_file)
{
    host_t<IFieldMemSp> init_perturb_mode(idx_range_kinsp);
    host_t<DFieldMemSp> init_perturb_amplitude(idx_range_kinsp);

    for (IdxSp const isp : idx_range_kinsp) {
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
        DFieldX const perturbation,
        int const perturb_mode,
        double const perturb_amplitude) const
{
    IdxRangeX const gridx = get_idx_range(perturbation);
    double const Lx = ddcHelper::total_interval_length(gridx);

    double const kx = perturb_mode * 2. * M_PI / Lx;
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            gridx,
            KOKKOS_LAMBDA(IdxX const ix) {
                CoordX const x = ddc::coordinate(ix);
                perturbation(ix) = perturb_amplitude * Kokkos::cos(kx * x);
            });
}
