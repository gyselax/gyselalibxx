// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include "chargedensitycalculator.hpp"

ChargeDensityCalculator::ChargeDensityCalculator(DViewVx coeffs) : m_quadrature(coeffs) {}

DSpanX ChargeDensityCalculator::operator()(DSpanX const rho, DViewSpXVx const allfdistribu) const
{
    Kokkos::Profiling::pushRegion("ChargeDensityCalculator");

    host_t<ViewSp<double>> const charges_host = ddc::host_discrete_space<IDimSp>().charges();
    IDomainSp const kin_species_domain = allfdistribu.domain<IDimSp>();
    host_t<ViewSp<double>> const kinetic_charges_host = charges_host[kin_species_domain];

    auto kinetic_charges_alloc
            = create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), kinetic_charges_host);
    ddc::ChunkSpan kinetic_charges = kinetic_charges_alloc.span_view();

    m_quadrature(
            Kokkos::DefaultExecutionSpace(),
            rho,
            KOKKOS_LAMBDA(IndexXVx ixvx) {
                double sum = 0.0;
                for (auto isp : kinetic_charges.domain()) {
                    sum += kinetic_charges(isp) * allfdistribu(isp, ixvx);
                }
                return sum;
            });

    IndexSp const last_kin_species = kin_species_domain.back();
    IndexSp const last_species = charges_host.domain().back();
    if (last_kin_species != last_species) {
        double const chargedens_adiabspecies = charge(last_species);

        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                rho.domain(),
                KOKKOS_LAMBDA(IndexX ix) { rho(ix) += chargedens_adiabspecies; });
    }

    Kokkos::Profiling::popRegion();

    return rho;
}
