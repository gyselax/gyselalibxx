// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include <ddc_helper.hpp>
#include <species_info.hpp>

#include "chargedensitycalculator.hpp"

ChargeDensityCalculator::ChargeDensityCalculator(const ChunkViewType& coeffs)
    : m_coefficients(coeffs)
{
}

void ChargeDensityCalculator::operator()(DSpanXY rho, DViewSpXYVxVy allfdistribu) const
{
    Kokkos::Profiling::pushRegion("ChargeDensityCalculator");
    IndexSp const last_kin_species = allfdistribu.domain<IDimSp>().back();
    IndexSp const last_species = ddc::discrete_space<IDimSp>().charges().domain().back();
    double chargedens_adiabspecies = 0.;
    if (last_kin_species != last_species) {
        chargedens_adiabspecies = double(charge(last_species));
    }

    // reduction over species and velocity space
    Kokkos::View<const double*****, Kokkos::LayoutRight> const allfdistribu_view
            = allfdistribu.allocation_kokkos_view();
    Kokkos::View<double**, Kokkos::LayoutRight> const rho_view = rho.allocation_kokkos_view();

    host_t<ViewSp<int>> const charges_host = ddc::host_discrete_space<IDimSp>().charges();
    host_t<ViewSp<int>> const kinetic_charges_host = charges_host[allfdistribu.domain<IDimSp>()];

    auto charges_alloc
            = create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), kinetic_charges_host);

    Kokkos::View<const int*, Kokkos::LayoutRight> const charges
            = charges_alloc.span_cview().allocation_kokkos_view();

    Kokkos::View<const double**, Kokkos::LayoutRight> const coef_view
            = m_coefficients.allocation_kokkos_view();

    std::size_t const nsp = allfdistribu_view.extent(0);
    std::size_t const nx = allfdistribu_view.extent(1);
    std::size_t const ny = allfdistribu_view.extent(2);
    std::size_t const nvx = allfdistribu_view.extent(3);
    std::size_t const nvy = allfdistribu_view.extent(4);

    Kokkos::parallel_for(
            Kokkos::TeamPolicy<>(nx * ny, Kokkos::AUTO),
            KOKKOS_LAMBDA(const Kokkos::TeamPolicy<>::member_type& team) {
                const int idx = team.league_rank();
                const int ix = idx / ny;
                const int iy = idx % ny;
                double teamSum = 0;
                Kokkos::parallel_reduce(
                        Kokkos::TeamThreadRange(team, nsp * nvx * nvy),
                        [&](int thread_idx, double& sum) {
                            int isp = thread_idx / (nvx * nvy);
                            thread_idx -= isp * (nvx * nvy);
                            int ivx = thread_idx / nvy;
                            int ivy = thread_idx % nvy;
                            sum += static_cast<double>(charges(isp)) * coef_view(ivx, ivy)
                                   * allfdistribu_view(isp, ix, iy, ivx, ivy);
                        },
                        teamSum);
                /*
                 * The code below is cleaner but currently only works on kokkos's develop branch (as of 15/03/24)
                 * The associated issue is : https://github.com/kokkos/kokkos/issues/6530
                 * After the next release of Kokkos this code should be uncommented to replace the more
                 * verbose version above.
                 *
                Kokkos::parallel_reduce(
                        Kokkos::TeamThreadMDRange(team, nsp, nvx, nvy),
                        [&](int isp, int ivx, int ivy, double& sum) {
                            sum += static_cast<double>(charges(isp)) * coef_view(ivx, ivy)
                                   * allfdistribu_view(isp, ix, iy, ivx, ivy);
                        },
                        teamSum);
                */
                rho_view(ix, iy) = chargedens_adiabspecies + teamSum;
            });
    Kokkos::Profiling::popRegion();
}
