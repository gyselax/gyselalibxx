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

void MpiMomentsCalculator::operator()(DFieldX rho, DConstFieldSpVxVyVzX allfdistribu) const
{
    Kokkos::Profiling::pushRegion("MpiChargeCalculator1d3v");

    DFieldMemX rho_local_alloc(get_idx_range(rho));
    DFieldX rho_local = get_field(rho_local_alloc);

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

void MpiMomentsCalculator::operator()(DFieldX kinetic, DConstFieldSpVxVyVzX allfdistribu, char ki) const
{
    Kokkos::Profiling::pushRegion("MpiKineticCalculator1d3v");

    DFieldMemX kinetic_local_alloc(get_idx_range(kinetic));
    DFieldX kinetic_local = get_field(kinetic_local_alloc);

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


void MpiMomentsCalculator::operator()(DFieldX mean_current_x, DFieldX mean_current_y, DFieldX mean_current_z, DFieldX rho, 
                                        DConstFieldSpVxVyVzX allfdistribu) const
{
    Kokkos::Profiling::pushRegion("MpiMeanCurrentCalculator1d3v");

    DFieldMemX mean_current_x_local_alloc(get_idx_range(mean_current_x));
    DFieldX mean_current_x_local = get_field(mean_current_x_local_alloc);

    DFieldMemX mean_current_y_local_alloc(get_idx_range(mean_current_y));
    DFieldX mean_current_y_local = get_field(mean_current_y_local_alloc);

    DFieldMemX mean_current_z_local_alloc(get_idx_range(mean_current_z));
    DFieldX mean_current_z_local = get_field(mean_current_z_local_alloc);

    m_local_moments_calculator(mean_current_x_local, mean_current_y_local, mean_current_z_local, rho, allfdistribu);

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

    MPI_Allreduce(
            mean_current_z_local.data_handle(),
            mean_current_z.data_handle(),
            mean_current_z.size(),
            MPI_type_descriptor_t<double>,
            MPI_SUM,
            m_comm);

    Kokkos::Profiling::popRegion();
}


void MpiMomentsCalculator::operator()(DFieldX momentum_x, DFieldX momentum_y, DFieldX momentum_z, 
    DConstFieldSpVxVyVzX allfdistribu) const
{
    Kokkos::Profiling::pushRegion("MpiMeanCurrentCalculator1d3v");

    DFieldMemX momentum_x_local_alloc(get_idx_range(momentum_x));
    DFieldX momentum_x_local = get_field(momentum_x_local_alloc);

    DFieldMemX momentum_y_local_alloc(get_idx_range(momentum_y));
    DFieldX momentum_y_local = get_field(momentum_y_local_alloc);

    DFieldMemX momentum_z_local_alloc(get_idx_range(momentum_z));
    DFieldX momentum_z_local = get_field(momentum_z_local_alloc);

    m_local_moments_calculator(momentum_x_local, momentum_y_local,  momentum_z_local, allfdistribu);

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

    MPI_Allreduce(
        momentum_z_local.data_handle(),
        momentum_z.data_handle(),
        momentum_z.size(),
        MPI_type_descriptor_t<double>,
        MPI_SUM,
        m_comm);

    Kokkos::Profiling::popRegion();
}


void MpiMomentsCalculator::operator()(DFieldSpX rho, DConstFieldSpVxVyVzX allfdistribu) const
{
    Kokkos::Profiling::pushRegion("MpiChargeCalculator1d3v");
    
    IdxRangeSp const kin_species_idx_range = get_idx_range<Species>(allfdistribu);

    DFieldMemSpX rho_local_alloc(get_idx_range(rho));
    DFieldSpX rho_local = get_field(rho_local_alloc);

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


void MpiMomentsCalculator::operator()(DFieldSpX mean_current_x, DFieldSpX mean_current_y,  DFieldSpX mean_current_z, DFieldSpX rho, 
                                        DConstFieldSpVxVyVzX allfdistribu) const
{
    Kokkos::Profiling::pushRegion("MpiMeanCurrentCalculator");

    IdxRangeSp const kin_species_idx_range = get_idx_range<Species>(allfdistribu);

    DFieldMemSpX mean_current_x_local_alloc(get_idx_range(mean_current_x));
    DFieldSpX mean_current_x_local = get_field(mean_current_x_local_alloc);

    DFieldMemSpX mean_current_y_local_alloc(get_idx_range(mean_current_y));
    DFieldSpX mean_current_y_local = get_field(mean_current_y_local_alloc);

    DFieldMemSpX mean_current_z_local_alloc(get_idx_range(mean_current_z));
    DFieldSpX mean_current_z_local = get_field(mean_current_z_local_alloc);

    m_local_moments_calculator(mean_current_x_local, mean_current_y_local, mean_current_z_local, rho, allfdistribu);

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

        MPI_Allreduce(
            mean_current_z_local[isp].data_handle(),
            mean_current_z[isp].data_handle(),
            mean_current_z[isp].size(),
            MPI_type_descriptor_t<double>,
            MPI_SUM,
            m_comm);
    }

    Kokkos::Profiling::popRegion();
}
