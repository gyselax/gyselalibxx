// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include "ddc_helper.hpp"
#include "singlemodeperturbinitialization.hpp"


SingleModePerturbInitialization::SingleModePerturbInitialization(
        DConstFieldSpVxVy fequilibrium,
        host_t<IFieldSp> init_perturb_mode,
        host_t<DFieldMemSp> init_perturb_amplitude)
    : m_fequilibrium(fequilibrium)
    , m_init_perturb_mode(std::move(init_perturb_mode))
    , m_init_perturb_amplitude(std::move(init_perturb_amplitude))
{
}

DFieldSpXYVxVy SingleModePerturbInitialization::operator()(DFieldSpXYVxVy const allfdistribu) const
{
    IdxRangeSp const gridsp = get_idx_range<Species>(allfdistribu);
    IdxRangeXY const gridxy = get_idx_range<GridX, GridY>(allfdistribu);
    IdxRangeXYVxVy const gridxyvxvy = get_idx_range<GridX, GridY, GridVx, GridVy>(allfdistribu);

    // Initialization of the perturbation
    DFieldMemXY perturbation_alloc(gridxy);
    ddc::ChunkSpan fequilibrium_proxy = get_field(m_fequilibrium);
    ddc::ChunkSpan perturbation_proxy = get_field(perturbation_alloc);
    ddc::for_each(gridsp, [&](IdxSp const isp) {
        perturbation_initialization(
                perturbation_proxy,
                m_init_perturb_mode(isp),
                m_init_perturb_amplitude(isp));

        // Initialization of the distribution function --> fill values
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                gridxyvxvy,
                KOKKOS_LAMBDA(IdxXYVxVy const ixyvxvy) {
                    IdxX const ix = ddc::select<GridX>(ixyvxvy);
                    IdxY const iy = ddc::select<GridY>(ixyvxvy);
                    IdxVx const ivx = ddc::select<GridVx>(ixyvxvy);
                    IdxVy const ivy = ddc::select<GridVy>(ixyvxvy);
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
        DConstFieldSpVxVy allfequilibrium,
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
        DFieldXY const perturbation,
        int const perturb_mode,
        double const perturb_amplitude) const
{
    IdxRangeXY const gridxy = get_idx_range<GridX, GridY>(perturbation);
    double const kx = perturb_mode * 2. * M_PI
                      / ddcHelper::total_interval_length(ddc::select<GridX>(gridxy));
    double const ky = perturb_mode * 2. * M_PI
                      / ddcHelper::total_interval_length(ddc::select<GridY>(gridxy));

    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            gridxy,
            KOKKOS_LAMBDA(IdxXY const ixy) {
                IdxX const ix = ddc::select<GridX>(ixy);
                CoordX const x = ddc::coordinate(ix);
                IdxY const iy = ddc::select<GridY>(ixy);
                CoordY const y = ddc::coordinate(iy);
                perturbation(ix, iy)
                        = perturb_amplitude * (Kokkos::cos(kx * x) + Kokkos::cos(ky * y));
            });
}
