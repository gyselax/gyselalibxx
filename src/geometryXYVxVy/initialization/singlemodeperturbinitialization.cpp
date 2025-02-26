// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_helper.hpp"
#include "singlemodeperturbinitialization.hpp"


SingleModePerturbInitialisation::SingleModePerturbInitialisation(
        DConstFieldSpVxVy fequilibrium,
        host_t<IFieldMemSp> init_perturb_mode,
        host_t<DFieldMemSp> init_perturb_amplitude)
    : m_fequilibrium(fequilibrium)
    , m_init_perturb_mode(std::move(init_perturb_mode))
    , m_init_perturb_amplitude(std::move(init_perturb_amplitude))
{
}

DFieldSpXYVxVy SingleModePerturbInitialisation::operator()(DFieldSpXYVxVy const allfdistribu) const
{
    IdxRangeSp const gridsp = get_idx_range<Species>(allfdistribu);
    IdxRangeXY const gridxy = get_idx_range<GridX, GridY>(allfdistribu);
    IdxRangeXYVxVy const gridxyvxvy = get_idx_range<GridX, GridY, GridVx, GridVy>(allfdistribu);

    // Initialisation of the perturbation
    DFieldMemXY perturbation_alloc(gridxy);
    DConstFieldSpVxVy fequilibrium_proxy = get_const_field(m_fequilibrium);
    DFieldXY perturbation_proxy = get_field(perturbation_alloc);
    ddc::for_each(gridsp, [&](IdxSp const isp) {
        perturbation_initialisation(
                perturbation_proxy,
                m_init_perturb_mode(isp),
                m_init_perturb_amplitude(isp));

        // Initialisation of the distribution function --> fill values
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


SingleModePerturbInitialisation SingleModePerturbInitialisation::init_from_input(
        DConstFieldSpVxVy allfequilibrium,
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

    return SingleModePerturbInitialisation(
            allfequilibrium,
            std::move(init_perturb_mode),
            std::move(init_perturb_amplitude));
}


void SingleModePerturbInitialisation::perturbation_initialisation(
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
