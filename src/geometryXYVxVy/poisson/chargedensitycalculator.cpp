// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include <ddc_helper.hpp>
#include <species_info.hpp>

#include "chargedensitycalculator.hpp"

ChargeDensityCalculator::ChargeDensityCalculator(DViewVxVy coeffs) : m_quadrature(coeffs) {}

void ChargeDensityCalculator::operator()(DSpanXY rho, DViewSpXYVxVy allfdistribu) const
{
    Kokkos::Profiling::pushRegion("ChargeDensityCalculator");

    IDomainSp const kin_species_domain = allfdistribu.domain<IDimSp>();
    host_t<DViewSp> const charges_host = ddc::host_discrete_space<IDimSp>().charges();
    host_t<DViewSp> const kinetic_charges_host = charges_host[allfdistribu.domain<IDimSp>()];

    auto const kinetic_charges_alloc
            = create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), kinetic_charges_host);
    ddc::ChunkSpan kinetic_charges = kinetic_charges_alloc.span_view();
    m_quadrature(
            Kokkos::DefaultExecutionSpace(),
            rho,
            KOKKOS_LAMBDA(IndexXYVxVy idx) {
                double sum = 0.0;
                for (auto isp : kinetic_charges.domain()) {
                    sum += kinetic_charges(isp) * allfdistribu(isp, idx);
                }
                return sum;
            });

    IndexSp const last_kin_species = kin_species_domain.back();
    IndexSp const last_species = charges_host.domain().back();
    if (last_kin_species != last_species) {
        double chargedens_adiabspecies = double(charge(last_species));
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                rho.domain(),
                KOKKOS_LAMBDA(IndexXY ixy) { rho(ixy) += chargedens_adiabspecies; });
    }

    Kokkos::Profiling::popRegion();
}
