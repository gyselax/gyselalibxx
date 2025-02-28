// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_helper.hpp"
#include "mpichargedensitycalculator.hpp"
#include "mpitools.hpp"
#include "species_info.hpp"

MpiChargeDensityCalculator::MpiChargeDensityCalculator(
        MPI_Comm comm,
        IChargeDensityCalculator const& local_charge_density_calculator)
    : m_local_charge_density_calculator(local_charge_density_calculator)
    , m_comm(comm)
{
}

void MpiChargeDensityCalculator::operator()(DFieldXY rho, DConstFieldSpXYVxVy allfdistribu) const
{
    Kokkos::Profiling::pushRegion("MpiChargeDensityCalculator");

    DFieldMemXY rho_local_alloc(get_idx_range(rho));
    DFieldXY rho_local = get_field(rho_local_alloc);

    m_local_charge_density_calculator(rho_local, allfdistribu);

    MPI_Allreduce(
            rho_local.data_handle(),
            rho.data_handle(),
            rho.size(),
            MPI_type_descriptor_t<double>,
            MPI_SUM,
            m_comm);

    Kokkos::Profiling::popRegion();
}
