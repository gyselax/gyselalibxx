// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_helper.hpp"
#include "mpimomentscalculator.hpp"
#include "mpitools.hpp"
#include "species_info.hpp"

MpiMomentsCalculator::MpiMomentsCalculator(
        MPI_Comm comm,
        IMomentsCalculator const& local_moments_calculator)
    : m_local_moments_calculator(local_moments_calculator)
    , m_comm(comm)
{
}

void MpiMomentsCalculator::operator()(DFieldXY rho, DConstFieldSpVxVyXY allfdistribu) const
{
    Kokkos::Profiling::pushRegion("MpiChargeCalculator");

    DFieldMemXY rho_local_alloc(get_idx_range(rho));
    DFieldXY rho_local = get_field(rho_local_alloc);

    m_local_moments_calculator(rho_local, allfdistribu);

    MPI_Allreduce(
            rho_local.data_handle(),
            rho.data_handle(),
            rho.size(),
            MPI_type_descriptor_t<double>,
            MPI_SUM,
            m_comm);

    Kokkos::Profiling::popRegion();
}

void MpiMomentsCalculator::operator()(DFieldXY kinetic, DConstFieldSpVxVyXY allfdistribu, char ki) const
{
    Kokkos::Profiling::pushRegion("MpiKineticCalculator");

    DFieldMemXY kinetic_local_alloc(get_idx_range(kinetic));
    DFieldXY kinetic_local = get_field(kinetic_local_alloc);

    m_local_moments_calculator(kinetic_local, allfdistribu, ki);

    MPI_Allreduce(
            kinetic_local.data_handle(),
            kinetic.data_handle(),
            kinetic.size(),
            MPI_type_descriptor_t<double>,
            MPI_SUM,
            m_comm);

    Kokkos::Profiling::popRegion();
}


void MpiMomentsCalculator::operator()(DFieldXY mean_current_x, DFieldXY mean_current_y, DFieldXY rho, 
                                        DConstFieldSpVxVyXY allfdistribu) const
{
    Kokkos::Profiling::pushRegion("MpiMeanCurrentCalculator");

    DFieldMemXY mean_current_x_local_alloc(get_idx_range(mean_current_x));
    DFieldXY mean_current_x_local = get_field(mean_current_x_local_alloc);

    DFieldMemXY mean_current_y_local_alloc(get_idx_range(mean_current_y));
    DFieldXY mean_current_y_local = get_field(mean_current_y_local_alloc);

    m_local_moments_calculator(mean_current_x_local, mean_current_y_local, rho, allfdistribu);

    MPI_Allreduce(
            mean_current_x_local.data_handle(),
            mean_current_x.data_handle(),
            mean_current_x.size(),
            MPI_type_descriptor_t<double>,
            MPI_SUM,
            m_comm);

    MPI_Allreduce(
            mean_current_y_local.data_handle(),
            mean_current_y.data_handle(),
            mean_current_y.size(),
            MPI_type_descriptor_t<double>,
            MPI_SUM,
            m_comm);

    Kokkos::Profiling::popRegion();
}


void MpiMomentsCalculator::operator()(DFieldXY momentum_x, DFieldXY momentum_y, 
    DConstFieldSpVxVyXY allfdistribu) const
{
    Kokkos::Profiling::pushRegion("MpiMeanCurrentCalculator");

    DFieldMemXY momentum_x_local_alloc(get_idx_range(momentum_x));
    DFieldXY momentum_x_local = get_field(momentum_x_local_alloc);

    DFieldMemXY momentum_y_local_alloc(get_idx_range(momentum_y));
    DFieldXY momentum_y_local = get_field(momentum_y_local_alloc);

    m_local_moments_calculator(momentum_x_local, momentum_y_local, allfdistribu);

    MPI_Allreduce(
        momentum_x_local.data_handle(),
        momentum_x.data_handle(),
        momentum_x.size(),
        MPI_type_descriptor_t<double>,
        MPI_SUM,
        m_comm);

    MPI_Allreduce(
        momentum_y_local.data_handle(),
        momentum_y.data_handle(),
        momentum_y.size(),
        MPI_type_descriptor_t<double>,
        MPI_SUM,
        m_comm);

    Kokkos::Profiling::popRegion();
}


void MpiMomentsCalculator::operator()(DFieldSpXY rho, DConstFieldSpVxVyXY allfdistribu) const
{
    Kokkos::Profiling::pushRegion("MpiChargeCalculator");
    
    IdxRangeSp const kin_species_idx_range = get_idx_range<Species>(allfdistribu);

    DFieldMemSpXY rho_local_alloc(get_idx_range(rho));
    DFieldSpXY rho_local = get_field(rho_local_alloc);

    m_local_moments_calculator(rho_local, allfdistribu);

    for (IdxSp isp : kin_species_idx_range) {
        MPI_Allreduce(
                rho_local[isp].data_handle(),
                rho[isp].data_handle(),
                rho[isp].size(),
                MPI_type_descriptor_t<double>,
                MPI_SUM,
                m_comm);
    }

    Kokkos::Profiling::popRegion();
}


void MpiMomentsCalculator::operator()(DFieldSpXY mean_current_x, DFieldSpXY mean_current_y, DFieldSpXY rho, 
                                        DConstFieldSpVxVyXY allfdistribu) const
{
    Kokkos::Profiling::pushRegion("MpiMeanCurrentCalculator");

    IdxRangeSp const kin_species_idx_range = get_idx_range<Species>(allfdistribu);

    DFieldMemSpXY mean_current_x_local_alloc(get_idx_range(mean_current_x));
    DFieldSpXY mean_current_x_local = get_field(mean_current_x_local_alloc);

    DFieldMemSpXY mean_current_y_local_alloc(get_idx_range(mean_current_y));
    DFieldSpXY mean_current_y_local = get_field(mean_current_y_local_alloc);

    m_local_moments_calculator(mean_current_x_local, mean_current_y_local, rho, allfdistribu);

    for (IdxSp isp : kin_species_idx_range) {
        MPI_Allreduce(
            mean_current_x_local[isp].data_handle(),
            mean_current_x[isp].data_handle(),
            mean_current_x[isp].size(),
            MPI_type_descriptor_t<double>,
            MPI_SUM,
            m_comm);

        MPI_Allreduce(
            mean_current_y_local[isp].data_handle(),
            mean_current_y[isp].data_handle(),
            mean_current_y[isp].size(),
            MPI_type_descriptor_t<double>,
            MPI_SUM,
            m_comm);
    }

    Kokkos::Profiling::popRegion();
}
