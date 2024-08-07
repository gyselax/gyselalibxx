// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include "chargedensitycalculator.hpp"

ChargeDensityCalculator::ChargeDensityCalculator(DConstFieldVx coeffs) : m_quadrature(coeffs) {}

DFieldX ChargeDensityCalculator::operator()(DFieldX const rho, DConstFieldSpXVx const allfdistribu)
        const
{
    Kokkos::Profiling::pushRegion("ChargeDensityCalculator");

    host_t<DConstFieldSp> const charges_host = ddc::host_discrete_space<Species>().charges();
    IdxRangeSp const kin_species_idx_range = get_idx_range<Species>(allfdistribu);
    host_t<DConstFieldSp> const kinetic_charges_host = charges_host[kin_species_idx_range];

    auto kinetic_charges_alloc
            = create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), kinetic_charges_host);
    ddc::ChunkSpan kinetic_charges = get_field(kinetic_charges_alloc);

    m_quadrature(
            Kokkos::DefaultExecutionSpace(),
            rho,
            KOKKOS_LAMBDA(IdxXVx ixvx) {
                double sum = 0.0;
                for (auto isp : get_idx_range(kinetic_charges)) {
                    sum += kinetic_charges(isp) * allfdistribu(isp, ixvx);
                }
                return sum;
            });

    IdxSp const last_kin_species = kin_species_idx_range.back();
    IdxSp const last_species = get_idx_range(charges_host).back();
    if (last_kin_species != last_species) {
        double const chargedens_adiabspecies = charge(last_species);

        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                get_idx_range(rho),
                KOKKOS_LAMBDA(IdxX ix) { rho(ix) += chargedens_adiabspecies; });
    }

    Kokkos::Profiling::popRegion();

    return rho;
}
