// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_helper.hpp"
#include "mpichargedensitycalculator.hpp"
#include "species_info.hpp"

MpiChargeDensityCalculator::MpiChargeDensityCalculator(
        IChargeDensityCalculator const& local_charge_density_calculator)
    : m_local_charge_density_calculator(local_charge_density_calculator)
{
}

void MpiChargeDensityCalculator::operator()(DFieldXY rho, DConstFieldSpVxVyXY allfdistribu) const
{
    Kokkos::Profiling::pushRegion("MpiChargeDensityCalculator");

    DFieldMemXY rho_local_alloc(
            "rho_local (MpiChargeDensityCalculator::operator())",
            get_idx_range(rho));
    DFieldXY rho_local = get_field(rho_local_alloc);

    m_local_charge_density_calculator(rho_local, allfdistribu);

    Kokkos::DefaultExecutionSpace().fence("Fence local ChargeDensityCalculator");

    ddc::parallel_deepcopy(rho, rho_local);

    Kokkos::Profiling::popRegion();
}
