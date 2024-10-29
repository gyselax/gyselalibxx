// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include "maxwellianequilibrium.hpp"

MaxwellianEquilibrium::MaxwellianEquilibrium(
        host_t<DFieldSpTor2D> density_eq_host,
        host_t<DFieldSpTor2D> temperature_eq_host,
        host_t<DFieldSpTor2D> mean_velocity_eq_host,
        host_t<DFieldTor2D> magnetic_field_host)
    : m_density_eq(get_idx_range(density_eq_host))
    , m_temperature_eq(get_idx_range(temperature_eq_host))
    , m_mean_velocity_eq(get_idx_range(mean_velocity_eq_host))
    , m_magnetic_field(get_idx_range(magnetic_field_host))
{
    ddc::parallel_deepcopy(m_density_eq, density_eq_host);
    ddc::parallel_deepcopy(m_temperature_eq, temperature_eq_host);
    ddc::parallel_deepcopy(m_mean_velocity_eq, mean_velocity_eq_host);
    ddc::parallel_deepcopy(m_magnetic_field, magnetic_field_host);
}

DFieldSpV2DTor2D MaxwellianEquilibrium::operator()(DFieldSpV2DTor2D const allfequilibrium) const
{
    IdxRangeSp const idxrange_sp(get_idx_range(allfequilibrium));
    IdxRangeV2DTor2D const idxrange_v2dtor2d(get_idx_range(allfequilibrium));

    // Initialization of the maxwellian
    ddc::for_each(idxrange_sp, [&](IdxSp const isp) {
        compute_maxwellian(
                allfequilibrium[isp],
                m_density_eq[isp],
                m_temperature_eq[isp],
                m_mean_velocity_eq[isp],
                m_magnetic_field);
    });
    return allfequilibrium;
}

void MaxwellianEquilibrium::compute_maxwellian(
        DFieldV2DTor2D const fMaxwellian,
        DConstFieldTor2D density,
        DConstFieldTor2D temperature,
        DConstFieldTor2D mean_velocity,
        DConstFieldTor2D magnetic_field)
{
    IdxRangeV2DTor2D const idxrange_v2dtor2d(get_idx_range(fMaxwellian));
    IdxRangeTor2D const idxrange_tor2d(get_idx_range(fMaxwellian));

    // Compute fM(vpar,mu) = (2*PI*T)**(-1.5)*n*exp(-energy) with
    //  - n the density, T the temperature and Upar the mean velocity
    //  - B the magnetic field and
    //  - energy = (0.5*(vpar-Upar)**2+mu*B)/T
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            idxrange_v2dtor2d,
            KOKKOS_LAMBDA(IdxV2DTor2D const iv2dtor2d) {
                IdxTor2D const idx_tor2d(iv2dtor2d);
                double const density_loc = density(idx_tor2d);
                double const temperature_loc = temperature(idx_tor2d);
                double const mean_velocity_loc = mean_velocity(idx_tor2d);
                double const magnetic_field_loc = magnetic_field(idx_tor2d);
                double const vpar = ddc::coordinate(ddc::select<GridVpar>(iv2dtor2d));
                double const mu = ddc::coordinate(ddc::select<GridMu>(iv2dtor2d));
                double const inv_2piT = 1. / (2. * M_PI * temperature_loc);
                double const coeff_maxw = Kokkos::sqrt(inv_2piT) * inv_2piT;
                double const energy = (0.5 * (vpar - mean_velocity_loc) * (vpar - mean_velocity_loc)
                                       + mu * magnetic_field_loc)
                                      / temperature_loc;
                fMaxwellian(iv2dtor2d) = coeff_maxw * density_loc * Kokkos::exp(-energy);
            });
}
