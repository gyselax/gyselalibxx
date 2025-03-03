// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include "chargedensitycalculator.hpp"
#include "ddc_alias_inline_functions.hpp"
#include "ddc_helper.hpp"
#include "species_info.hpp"

ChargeDensityCalculator::ChargeDensityCalculator(DConstFieldVxVy coeffs) : m_quadrature(coeffs) {}

void ChargeDensityCalculator::operator()(DFieldXY rho, DConstFieldSpVxVyXY allfdistribu) const
{
    Kokkos::Profiling::pushRegion("ChargeDensityCalculator");

    IdxRangeSp const kin_species_idx_range = get_idx_range<Species>(allfdistribu);
    host_t<DConstFieldSp> const charges_host = ddc::host_discrete_space<Species>().charges();
    host_t<DConstFieldSp> const kinetic_charges_host
            = charges_host[get_idx_range<Species>(allfdistribu)];

    auto const kinetic_charges_alloc
            = create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), kinetic_charges_host);
    DConstFieldSp kinetic_charges = get_const_field(kinetic_charges_alloc);
    m_quadrature(
            Kokkos::DefaultExecutionSpace(),
            rho,
            KOKKOS_LAMBDA(IdxXYVxVy idx) {
                double sum = 0.0;
                for (auto isp : kin_species_idx_range) {
                    sum += kinetic_charges(isp) * allfdistribu(isp, idx);
                }
                return sum;
            });

    IdxSp const last_kin_species = kin_species_idx_range.back();
    IdxSp const last_species = get_idx_range(charges_host).back();
    if (last_kin_species != last_species) {
        double chargedens_adiabspecies = double(charge(last_species));
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                get_idx_range(rho),
                KOKKOS_LAMBDA(IdxXY ixy) { rho(ixy) += chargedens_adiabspecies; });
    }

    Kokkos::Profiling::popRegion();
}
