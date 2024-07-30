// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include "bumpontailequilibrium.hpp"

BumpontailEquilibrium::BumpontailEquilibrium(
        host_t<DFieldSp> epsilon_bot,
        host_t<DFieldSp> temperature_bot,
        host_t<DFieldSp> mean_velocity_bot)
    : m_epsilon_bot(std::move(epsilon_bot))
    , m_temperature_bot(std::move(temperature_bot))
    , m_mean_velocity_bot(std::move(mean_velocity_bot))
{
}

DSpanSpVx BumpontailEquilibrium::operator()(DSpanSpVx const allfequilibrium) const
{
    IDomainVx const gridvx = allfequilibrium.domain<IDimVx>();
    IDomainSp const gridsp = allfequilibrium.domain<IDimSp>();

    // Initialization of the maxwellian
    DFieldVx maxwellian_alloc(gridvx);
    ddc::ChunkSpan maxwellian = maxwellian_alloc.span_view();
    ddc::for_each(gridsp, [&](IndexSp const isp) {
        compute_twomaxwellian(
                maxwellian,
                m_epsilon_bot(isp),
                m_temperature_bot(isp),
                m_mean_velocity_bot(isp));

        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                gridvx,
                KOKKOS_LAMBDA(IndexVx const ivx) { allfequilibrium(isp, ivx) = maxwellian(ivx); });
    });
    return allfequilibrium;
}

BumpontailEquilibrium BumpontailEquilibrium::init_from_input(
        IDomainSp dom_kinsp,
        PC_tree_t const& yaml_input_file)
{
    host_t<DFieldSp> epsilon_bot(dom_kinsp);
    host_t<DFieldSp> temperature_bot(dom_kinsp);
    host_t<DFieldSp> mean_velocity_bot(dom_kinsp);

    for (IndexSp const isp : dom_kinsp) {
        PC_tree_t const conf_isp = PCpp_get(yaml_input_file, ".SpeciesInfo[%d]", isp.uid());

        epsilon_bot(isp) = PCpp_double(conf_isp, ".epsilon_bot");
        temperature_bot(isp) = PCpp_double(conf_isp, ".temperature_bot");
        mean_velocity_bot(isp) = PCpp_double(conf_isp, ".mean_velocity_bot");
    }

    return BumpontailEquilibrium(
            std::move(epsilon_bot),
            std::move(temperature_bot),
            std::move(mean_velocity_bot));
}


void BumpontailEquilibrium::compute_twomaxwellian(
        DSpanVx const fMaxwellian,
        double const epsilon_bot,
        double const temperature_bot,
        double const mean_velocity_bot) const
{
    double const inv_sqrt_2pi = 1. / sqrt(2. * M_PI);
    double const norm_f2 = inv_sqrt_2pi / sqrt(temperature_bot);
    IDomainVx const gridvx = fMaxwellian.domain();
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            gridvx,
            KOKKOS_LAMBDA(IndexVx const ivx) {
                CoordVx const vx = ddc::coordinate(ivx);
                // bulk plasma particles
                double const f1_v = (1. - epsilon_bot) * inv_sqrt_2pi * Kokkos::exp(-0.5 * vx * vx);
                // beam
                double const f2_v = epsilon_bot * norm_f2
                                    * (Kokkos::exp(
                                            -(vx - mean_velocity_bot) * (vx - mean_velocity_bot)
                                            / (2. * temperature_bot)));
                // fM(v) = f1(v) + f2(v)
                fMaxwellian(ivx) = f1_v + f2_v;
            });
}
