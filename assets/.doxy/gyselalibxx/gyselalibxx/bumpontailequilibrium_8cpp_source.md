

# File bumpontailequilibrium.cpp

[**File List**](files.md) **>** [**geometryXVx**](dir_e51b496b46dd687775e46e0826614574.md) **>** [**initialisation**](dir_cdb336346544d0d5f695f9cdfe73a70e.md) **>** [**bumpontailequilibrium.cpp**](bumpontailequilibrium_8cpp.md)

[Go to the documentation of this file](bumpontailequilibrium_8cpp.md)


```C++
// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include "bumpontailequilibrium.hpp"

BumpontailEquilibrium::BumpontailEquilibrium(
        host_t<DFieldMemSp> epsilon_bot,
        host_t<DFieldMemSp> temperature_bot,
        host_t<DFieldMemSp> mean_velocity_bot)
    : m_epsilon_bot(std::move(epsilon_bot))
    , m_temperature_bot(std::move(temperature_bot))
    , m_mean_velocity_bot(std::move(mean_velocity_bot))
{
}

DFieldSpVx BumpontailEquilibrium::operator()(DFieldSpVx const allfequilibrium) const
{
    IdxRangeVx const gridvx = get_idx_range<GridVx>(allfequilibrium);
    IdxRangeSp const gridsp = get_idx_range<Species>(allfequilibrium);

    // Initialisation of the maxwellian
    DFieldMemVx maxwellian_alloc(gridvx);
    DFieldVx maxwellian = get_field(maxwellian_alloc);
    ddc::for_each(gridsp, [&](IdxSp const isp) {
        compute_twomaxwellian(
                maxwellian,
                m_epsilon_bot(isp),
                m_temperature_bot(isp),
                m_mean_velocity_bot(isp));

        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                gridvx,
                KOKKOS_LAMBDA(IdxVx const ivx) { allfequilibrium(isp, ivx) = maxwellian(ivx); });
    });
    return allfequilibrium;
}

BumpontailEquilibrium BumpontailEquilibrium::init_from_input(
        IdxRangeSp idx_range_kinsp,
        PC_tree_t const& yaml_input_file)
{
    host_t<DFieldMemSp> epsilon_bot(idx_range_kinsp);
    host_t<DFieldMemSp> temperature_bot(idx_range_kinsp);
    host_t<DFieldMemSp> mean_velocity_bot(idx_range_kinsp);

    for (IdxSp const isp : idx_range_kinsp) {
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
        DFieldVx const fMaxwellian,
        double const epsilon_bot,
        double const temperature_bot,
        double const mean_velocity_bot) const
{
    double const inv_sqrt_2pi = 1. / sqrt(2. * M_PI);
    double const norm_f2 = inv_sqrt_2pi / sqrt(temperature_bot);
    IdxRangeVx const gridvx = get_idx_range(fMaxwellian);
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            gridvx,
            KOKKOS_LAMBDA(IdxVx const ivx) {
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
```


