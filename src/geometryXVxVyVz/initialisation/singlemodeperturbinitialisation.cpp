// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_helper.hpp"
#include "singlemodeperturbinitialisation.hpp"


SingleModePerturbInitialisation::SingleModePerturbInitialisation(
        DConstFieldSpVxVyVz fequilibrium,
        host_t<IFieldMemSp> init_perturb_mode,
        host_t<DFieldMemSp> init_perturb_amplitude)
    : m_fequilibrium(fequilibrium)
    , m_init_perturb_mode(std::move(init_perturb_mode))
    , m_init_perturb_amplitude(std::move(init_perturb_amplitude))
{
}

DFieldSpXVxVyVz SingleModePerturbInitialisation::operator()(DFieldSpXVxVyVz const allfdistribu) const
{
    IdxRangeSp const gridsp = get_idx_range<Species>(allfdistribu);
    IdxRangeX const gridx = get_idx_range<GridX>(allfdistribu);
    IdxRangeXVxVyVz const gridxvxvyvz = get_idx_range<GridX, GridVx, GridVy, GridVz>(allfdistribu);

    // Initialisation of the perturbation
    DFieldMemX perturbation_alloc(gridx);
    DConstFieldSpVxVyVz fequilibrium_proxy = get_const_field(m_fequilibrium);
    DFieldX perturbation_proxy = get_field(perturbation_alloc);
    ddc::for_each(gridsp, [&](IdxSp const isp) {
        perturbation_initialisation(
                perturbation_proxy,
                m_init_perturb_mode(isp),
                m_init_perturb_amplitude(isp));

        double w0 = 0.17773094298487538;
        double k0 = 0.196;
        double temperature = 0.125;
        double const inv_2pi = 1. / (2. * M_PI * temperature);

        // Initialisation of the distribution function --> fill values
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                gridxvxvyvz,
                KOKKOS_LAMBDA(IdxXVxVyVz const ixvxvyvz) {
                    /*
                    IdxX const ix = ddc::select<GridX>(ixvxvyvz);
                    IdxVx const ivx = ddc::select<GridVx>(ixvxvyvz);
                    IdxVy const ivy = ddc::select<GridVy>(ixvxvyvz);
                    IdxVz const ivz = ddc::select<GridVz>(ixvxvyvz);
                    double fdistribu_val
                            = fequilibrium_proxy(isp, ivx, ivy, ivz) * (1. + perturbation_proxy(ix));
                    if (fdistribu_val < 1.e-60) {
                        fdistribu_val = 1.e-60;
                    }
                    allfdistribu(isp, ix, ivx, ivy, ivz) = fdistribu_val;
                    */
                    IdxX const ix = ddc::select<GridX>(ixvxvyvz);
                    IdxVx const ivx = ddc::select<GridVx>(ixvxvyvz);
                    IdxVy const ivy = ddc::select<GridVy>(ixvxvyvz);
                    IdxVz const ivz = ddc::select<GridVz>(ixvxvyvz);
                    double const x = ddc::coordinate(ddc::select<GridX>(ixvxvyvz));
                    double const vx = ddc::coordinate(ddc::select<GridVx>(ixvxvyvz));
                    double const vy = ddc::coordinate(ddc::select<GridVy>(ixvxvyvz));
                    double const vz = ddc::coordinate(ddc::select<GridVz>(ixvxvyvz));

                    double fdistribu_val = std::pow(inv_2pi, 1.5) 
                                     * Kokkos::exp(
                                             -(vx*vx)
                                             / (2.0*temperature))
                                     * Kokkos::exp(
                                             -((vy - w0/k0 * Kokkos::sin(0.196 * x)) * (vy - w0/k0 * Kokkos::sin(0.196 * x))
                                             + (vz + w0/k0 * Kokkos::cos(0.196 * x)) * (vz + w0/k0 * Kokkos::cos(0.196 * x)))
                                             / (2.0*temperature ));

                    if (fdistribu_val < 1.e-60) {
                        fdistribu_val = 1.e-60;
                    }
                    allfdistribu(isp, ix, ivx, ivy, ivz) = fdistribu_val;
                });
    });
    return allfdistribu;
}


SingleModePerturbInitialisation SingleModePerturbInitialisation::init_from_input(
        DConstFieldSpVxVyVz allfequilibrium,
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
        DFieldX const perturbation,
        int const perturb_mode,
        double const perturb_amplitude) const
{
    IdxRangeX const gridx = get_idx_range<GridX>(perturbation);
    double const kx = perturb_mode * 2. * M_PI
                      / ddcHelper::total_interval_length(ddc::select<GridX>(gridx));

    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            gridx,
            KOKKOS_LAMBDA(IdxX const ix) {
                IdxX const ixx = ddc::select<GridX>(ix);
                CoordX const x = ddc::coordinate(ix);
               
                perturbation(ix)
                        = (0.0 * perturb_amplitude) * Kokkos::sin(kx * x);
            });
}
