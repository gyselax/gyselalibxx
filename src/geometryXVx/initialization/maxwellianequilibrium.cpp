// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include "maxwellianequilibrium.hpp"

MaxwellianEquilibrium::MaxwellianEquilibrium(
        host_t<DFieldSp> density_eq,
        host_t<DFieldSp> temperature_eq,
        host_t<DFieldSp> mean_velocity_eq)
    : m_density_eq(std::move(density_eq))
    , m_temperature_eq(std::move(temperature_eq))
    , m_mean_velocity_eq(std::move(mean_velocity_eq))
{
}

DSpanSpVx MaxwellianEquilibrium::operator()(DSpanSpVx const allfequilibrium) const
{
    IDomainVx const gridvx = allfequilibrium.domain<IDimVx>();
    IDomainSp const gridsp = allfequilibrium.domain<IDimSp>();

    // Initialization of the maxwellian
    DFieldVx maxwellian_alloc(gridvx);
    ddc::ChunkSpan maxwellian = maxwellian_alloc.span_view();
    ddc::for_each(gridsp, [&](IndexSp const isp) {
        compute_maxwellian(
                maxwellian,
                m_density_eq(isp),
                m_temperature_eq(isp),
                m_mean_velocity_eq(isp));

        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                gridvx,
                KOKKOS_LAMBDA(IndexVx const ivx) { allfequilibrium(isp, ivx) = maxwellian(ivx); });
    });
    return allfequilibrium;
}

void MaxwellianEquilibrium::compute_maxwellian(
        DSpanVx const fMaxwellian,
        double const density,
        double const temperature,
        double const mean_velocity)
{
    double const inv_sqrt_2piT = 1. / Kokkos::sqrt(2. * M_PI * temperature);
    IDomainVx const gridvx = fMaxwellian.domain();
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            gridvx,
            KOKKOS_LAMBDA(IndexVx const ivx) {
                CoordVx const vx = ddc::coordinate(ivx);
                fMaxwellian(ivx) = density * inv_sqrt_2piT
                                   * Kokkos::exp(
                                           -(vx - mean_velocity) * (vx - mean_velocity)
                                           / (2. * temperature));
            });
}
