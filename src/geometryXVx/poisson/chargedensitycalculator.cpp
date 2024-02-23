// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include "chargedensitycalculator.hpp"

ChargeDensityCalculator::ChargeDensityCalculator(const ChunkViewType& coeffs)
    : m_coefficients(coeffs)
{
}

device_t<DSpanX> ChargeDensityCalculator::operator()(
        device_t<DSpanX> const rho,
        device_t<DViewSpXVx> const allfdistribu) const
{
    Kokkos::Profiling::pushRegion("ChargeDensityCalculator");
    IndexSp const last_kin_species = allfdistribu.domain<IDimSp>().back();
    IndexSp const last_species = ddc::discrete_space<IDimSp>().charges().domain().back();
    double chargedens_adiabspecies = 0.;
    if (last_kin_species != last_species) {
        chargedens_adiabspecies = double(charge(last_species));
    }

    // reduction over species and velocity space
    Kokkos::View<const double***, Kokkos::LayoutRight> const allfdistribu_view
            = allfdistribu.allocation_kokkos_view();
    Kokkos::View<double*, Kokkos::LayoutRight> const rho_view = rho.allocation_kokkos_view();

    ViewSp<int> const charges = ddc::host_discrete_space<IDimSp>().charges();
    ViewSp<int> const kinetic_charges = charges[allfdistribu.domain<IDimSp>()];

    auto charges_alloc
            = create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), kinetic_charges);

    Kokkos::View<const int*, Kokkos::LayoutRight> const charges_view
            = charges_alloc.span_cview().allocation_kokkos_view();

    Kokkos::View<const double*, Kokkos::LayoutRight> const coef_view
            = m_coefficients.allocation_kokkos_view();

    std::size_t const nsp = allfdistribu_view.extent(0);
    std::size_t const nx = allfdistribu_view.extent(1);
    std::size_t const nvx = allfdistribu_view.extent(2);

    Kokkos::parallel_for(
            Kokkos::TeamPolicy<>(nx, Kokkos::AUTO),
            KOKKOS_LAMBDA(const Kokkos::TeamPolicy<>::member_type& team) {
                const int ix = team.league_rank();
                double teamSum = 0;
                Kokkos::parallel_reduce(
                        Kokkos::TeamThreadRange(team, nvx),
                        [&](int const& ivx, double& sum) {
                            // [TO DO] Nested reduction may be possible?
                            for (std::size_t isp = 0; isp < nsp; isp++) {
                                sum += static_cast<double>(charges_view(isp)) * coef_view(ivx)
                                       * allfdistribu_view(isp, ix, ivx);
                            }
                        },
                        teamSum);
                rho_view(ix) = chargedens_adiabspecies + teamSum;
            });
    Kokkos::Profiling::popRegion();

    return rho;
}
