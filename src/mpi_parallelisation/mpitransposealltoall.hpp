// SPDX-License-Identifier: MIT
#pragma once
#include <numeric>

#include <ddc/ddc.hpp>

#include "ddc_helper.hpp"
#include "impitranspose.hpp"
#include "mpilayout.hpp"
#include "mpitools.hpp"
#include "transpose.hpp"

/**
 * @brief A class describing an operator for converting from/to different MPI layouts using AlltoAll.
 *
 * This class implements a basic AlltoAll operator and currently only works with a basic MPIBlockLayout.
 *
 * @tparam Layout1 One of the MPI layouts.
 * @tparam Layout2 The other MPI layouts.
 */
template <class Layout1, class Layout2>
class MPITransposeAllToAll : public IMPITranspose<Layout1, Layout2>
{
public:
    /// The type of the domain of the first MPI layout.
    using domain_type1 = typename Layout1::discrete_domain_type;
    /// The type of the domain of the second MPI layout.
    using domain_type2 = typename Layout2::discrete_domain_type;
    /// The type of the domain of the first MPI layout.
    using distributed_domain_type1 = typename Layout1::distributed_subdomain;
    /// The type of the domain of the second MPI layout.
    using distributed_domain_type2 = typename Layout2::distributed_subdomain;
    /// The way in which data should be laid out in memory on Chunks being transposed.
    using chunk_layout_type = std::experimental::layout_right;

private:
    using layout_1_mpi_dims = ddcHelper::
            apply_template_to_type_seq_t<MPIDim, typename Layout1::distributed_type_seq>;
    using layout_2_mpi_dims = ddcHelper::
            apply_template_to_type_seq_t<MPIDim, typename Layout2::distributed_type_seq>;
    using layout_1_mpi_domain_type
            = ddc::detail::convert_type_seq_to_discrete_domain<layout_1_mpi_dims>;
    using layout_2_mpi_domain_type
            = ddc::detail::convert_type_seq_to_discrete_domain<layout_2_mpi_dims>;

private:
    int m_comm_size;
    Layout1 m_layout_1;
    Layout2 m_layout_2;
    domain_type1 m_local_domain_1;
    domain_type2 m_local_domain_2;
    layout_1_mpi_domain_type m_layout_1_mpi_domain;
    layout_2_mpi_domain_type m_layout_2_mpi_domain;

public:
    /**
     * @brief A constructor for the transpose operator.
     *
     * @param global_domain The global domain of the data across processes provided in either layout.
     * @param comm The MPI communicator.
     */
    template <class Domain>
    MPITransposeAllToAll(Domain global_domain, MPI_Comm comm)
        : IMPITranspose<Layout1, Layout2>(comm)
    {
        static_assert(
                std::is_same_v<Domain, domain_type1> || std::is_same_v<Domain, domain_type2>,
                "The initialisation global domain should be described by one of the layouts");
        domain_type1 global_domain_layout_1(global_domain);
        domain_type2 global_domain_layout_2(global_domain);
        int rank;
        MPI_Comm_size(comm, &m_comm_size);
        MPI_Comm_rank(comm, &rank);
        distributed_domain_type1 distrib_domain(global_domain_layout_1);
        int n_cores_max = distrib_domain.size();
        assert(m_comm_size <= n_cores_max);
        int n_elems = n_cores_max / m_comm_size;
        // Ensure that the load balancing is good
        assert(n_cores_max - m_comm_size * n_elems == 0);
        m_local_domain_1 = m_layout_1.distribute_domain(global_domain_layout_1, m_comm_size, rank);
        m_local_domain_2 = m_layout_2.distribute_domain(global_domain_layout_2, m_comm_size, rank);
        m_layout_1_mpi_domain = get_distribution(
                distributed_domain_type1(m_local_domain_1),
                distributed_domain_type1(global_domain_layout_1));
        m_layout_2_mpi_domain = get_distribution(
                distributed_domain_type2(m_local_domain_2),
                distributed_domain_type2(global_domain_layout_2));
    }

    /**
     * @brief Getter for the local domain.
     *
     * @tparam Layout The layout whose domain should be retrieved.
     *
     * @returns The local domain in the specified MPI layout.
     */
    template <class Layout>
    auto get_local_domain()
    {
        static_assert(
                std::is_same_v<Layout, Layout1> || std::is_same_v<Layout, Layout2>,
                "Transpose class does not handle requested layout");
        if constexpr (std::is_same_v<Layout, Layout1>) {
            return m_local_domain_1;
        } else if constexpr (std::is_same_v<Layout, Layout2>) {
            return m_local_domain_2;
        }
    }

    /**
     * @brief An operator which transposes from one layout to another.
     *
     * If the dimensions are ordered differently in the domains of the two layouts
     * then this function can be used. If the domains have the same type then the
     * transpose_to function must be used to disambiguate.
     *
     * @param[in] execution_space The execution space (Host/Device) where the code will run.
     * @param[out] recv_span The chunk which will describe the data in the new layout. This
     *                      data is gathered from the MPI processes.
     * @param[in] send_span The chunk describing the data in the current layout. This data
     *                      will be scattered to other MPI processes.
     */
    template <class ElementType, class InDomain, class OutDomain, class MemSpace, class ExecSpace>
    void operator()(
            ExecSpace const& execution_space,
            ddc::ChunkSpan<ElementType, OutDomain, chunk_layout_type, MemSpace> recv_span,
            ddc::ChunkView<ElementType, InDomain, chunk_layout_type, MemSpace> send_span)
    {
        static_assert(!std::is_same_v<InDomain, OutDomain>);
        static_assert(
                (std::is_same_v<InDomain, typename Layout1::discrete_domain_type>)
                || (std::is_same_v<InDomain, typename Layout2::discrete_domain_type>));
        static_assert(
                (std::is_same_v<OutDomain, typename Layout1::discrete_domain_type>)
                || (std::is_same_v<OutDomain, typename Layout2::discrete_domain_type>));
        using InLayout = std::conditional_t<
                std::is_same_v<InDomain, typename Layout1::discrete_domain_type>,
                Layout1,
                Layout2>;
        using OutLayout = std::conditional_t<
                std::is_same_v<OutDomain, typename Layout1::discrete_domain_type>,
                Layout1,
                Layout2>;
        this->transpose_to<OutLayout>(execution_space, recv_span, send_span);
    }

    /**
     * @brief An operator which transposes from one layout to another.
     *
     * If the dimensions are ordered differently in the domains of the two layouts
     * then this function can be used. If the domains have the same type then the
     * transpose_to function must be used to disambiguate.
     *
     * @tparam OutLayout The layout that the data should be transposed to.
     *
     * @param[in] execution_space The execution space (Host/Device) where the code will run.
     * @param[out] recv_span The chunk which will describe the data in the new layout. This
     *                      data is gathered from the MPI processes.
     * @param[in] send_span The chunk describing the data in the current layout. This data
     *                      will be scattered to other MPI processes.
     */
    template <class OutLayout, class ElementType, class MemSpace, class ExecSpace, class InDomain>
    void transpose_to(
            ExecSpace const& execution_space,
            ddc::ChunkSpan<
                    ElementType,
                    typename OutLayout::discrete_domain_type,
                    chunk_layout_type,
                    MemSpace> recv_span,
            ddc::ChunkView<ElementType, InDomain, chunk_layout_type, MemSpace> send_span)
    {
        using InLayout = std::conditional_t<std::is_same_v<OutLayout, Layout1>, Layout2, Layout1>;
        /*****************************************************************
         * Validate input
         *****************************************************************/
        static_assert((std::is_same_v<InLayout, Layout1>) || (std::is_same_v<InLayout, Layout2>));
        static_assert((std::is_same_v<OutLayout, Layout1>) || (std::is_same_v<OutLayout, Layout2>));
        static_assert(Kokkos::SpaceAccessibility<ExecSpace, MemSpace>::accessible);

        static_assert(std::is_same_v<InDomain, typename InLayout::discrete_domain_type>);
        using OutDomain = typename OutLayout::discrete_domain_type;

        /*****************************************************************
         * Define groups of tags to build necessary transform
         *****************************************************************/
        // Currently distributed dims will be gathered
        using dims_to_gather = typename InLayout::distributed_type_seq;
        // Future distributed dims must be scattered
        using dims_to_scatter = typename OutLayout::distributed_type_seq;
        // Get MPI tags (e.g. MPI<DimX> indicates where DimX will be gathered from)
        using gather_mpi_dims = ddcHelper::
                apply_template_to_type_seq_t<MPIDim, typename InLayout::distributed_type_seq>;
        using scatter_mpi_dims = ddcHelper::
                apply_template_to_type_seq_t<MPIDim, typename OutLayout::distributed_type_seq>;
        // Get complete set of tags from domains
        using input_ordered_dims = ddc::to_type_seq_t<InDomain>;
        using output_ordered_dims = ddc::to_type_seq_t<OutDomain>;
        // Find tags that are neither scattered nor gathered
        using batch_dims = ddc::type_seq_remove_t<
                ddc::type_seq_remove_t<input_ordered_dims, dims_to_scatter>,
                dims_to_gather>;
        // Insert the MPI tags into the received data to introduce the dimension split
        using input_mpi_domain_tags
                = insert_mpi_tags_into_seq_t<scatter_mpi_dims, input_ordered_dims>;
        using output_mpi_domain_tags
                = insert_mpi_tags_into_seq_t<gather_mpi_dims, output_ordered_dims>;
        // During the alltoall call the MPI tags must be first to correctly send/receive data
        // but data does not reorder itself
        using input_alltoall_dim_order
                = ddc::type_seq_merge_t<scatter_mpi_dims, input_ordered_dims>;
        using output_alltoall_dim_order
                = ddc::type_seq_merge_t<gather_mpi_dims, input_ordered_dims>;

        /*****************************************************************
         * Define useful domain types
         *****************************************************************/
        // Get the domains of objects which are distributed across processes
        using gather_domain_type = typename InLayout::distributed_subdomain;
        // Get the domains of objects which will be distributed across processes
        using scatter_domain_type = typename OutLayout::distributed_subdomain;
        // Get the domain of objects that are not distributed across processes
        using batch_domain_type = ddc::detail::convert_type_seq_to_discrete_domain<batch_dims>;
        using gather_mpi_domain_type
                = ddc::detail::convert_type_seq_to_discrete_domain<gather_mpi_dims>;
        using scatter_mpi_domain_type
                = ddc::detail::convert_type_seq_to_discrete_domain<scatter_mpi_dims>;
        // Get the domain containing MPI tags which can be used to describe the function input
        using input_mpi_domain_type
                = ddc::detail::convert_type_seq_to_discrete_domain<input_mpi_domain_tags>;
        // Get the domain containing MPI tags which can be used to describe the function output
        using output_mpi_domain_type
                = ddc::detail::convert_type_seq_to_discrete_domain<output_mpi_domain_tags>;
        // Get the domain that should be used when scattering data
        using input_alltoall_domain_type
                = ddc::detail::convert_type_seq_to_discrete_domain<input_alltoall_dim_order>;
        // Get the domain that should be used when gathering data
        using output_alltoall_domain_type
                = ddc::detail::convert_type_seq_to_discrete_domain<output_alltoall_dim_order>;

        /*****************************************************************
         * Build domains
         *****************************************************************/
        // Collect the artificial dimension describing the MPI rank where the scattered information
        // will be sent to or where the gathered information will be collected from
        gather_mpi_domain_type gather_mpi_domain;
        scatter_mpi_domain_type scatter_mpi_domain;
        if constexpr (std::is_same_v<InLayout, Layout1>) {
            gather_mpi_domain = m_layout_1_mpi_domain;
            scatter_mpi_domain = m_layout_2_mpi_domain;
        } else {
            gather_mpi_domain = m_layout_2_mpi_domain;
            scatter_mpi_domain = m_layout_1_mpi_domain;
        }

        // Collect the useful subdomains described in the fields
        batch_domain_type batch_domain(send_span.domain());
        scatter_domain_type scatter_domain(recv_span.domain());
        gather_domain_type gather_domain(send_span.domain());

        // Build the domains describing the function inputs but including the MPI rank information
        input_mpi_domain_type
                input_mpi_domain(scatter_mpi_domain, scatter_domain, gather_domain, batch_domain);
        output_mpi_domain_type
                output_mpi_domain(gather_mpi_domain, scatter_domain, gather_domain, batch_domain);
        assert(input_mpi_domain.size() == send_span.size());
        assert(output_mpi_domain.size() == recv_span.size());

        // Create the domains used during the alltoall call (the MPI rank domain is first in this layout)
        input_alltoall_domain_type input_alltoall_domain(input_mpi_domain);
        output_alltoall_domain_type output_alltoall_domain(output_mpi_domain);

        /*****************************************************************
         * Create views on the function inputs with the domains including the MPI rank information
         *****************************************************************/
        ddc::ChunkView<ElementType, input_mpi_domain_type, chunk_layout_type, MemSpace>
                send_mpi_span(send_span.data_handle(), input_mpi_domain);
        ddc::ChunkSpan<ElementType, output_mpi_domain_type, chunk_layout_type, MemSpace>
                recv_mpi_span(recv_span.data_handle(), output_mpi_domain);

        /*****************************************************************
         * Transpose data (both on the rank and between ranks)
         *****************************************************************/
        // Create views or copies of the function inputs such that they are laid out on the
        // domain used during the alltoall call
        auto alltoall_send_buffer = ddcHelper::create_transpose_mirror_view_and_copy<
                input_alltoall_domain_type>(execution_space, send_mpi_span);
        auto alltoall_recv_buffer = ddcHelper::create_transpose_mirror<
                output_alltoall_domain_type>(execution_space, recv_mpi_span);

        // Call the MPI AlltoAll routine
        call_all_to_all(
                execution_space,
                alltoall_recv_buffer.span_view(),
                alltoall_send_buffer.span_cview());

        // If alltoall_recv_buffer owns its data (not a view) then copy the results back to
        // recv_mpi_span which is a view on recv_span, the function output
        if constexpr (!ddc::is_borrowed_chunk_v<decltype(alltoall_recv_buffer)>) {
            transpose_layout(execution_space, recv_mpi_span, alltoall_recv_buffer.span_cview());
        }
    }

private:
    /// Function handling the MPI call
    template <
            class ElementType,
            class MPIRecvDomain,
            class MPISendDomain,
            class MemSpace,
            class ExecSpace>
    void call_all_to_all(
            ExecSpace const& execution_space,
            ddc::ChunkSpan<ElementType, MPIRecvDomain, chunk_layout_type, MemSpace> recv_span,
            ddc::ChunkView<ElementType, MPISendDomain, chunk_layout_type, MemSpace> send_span)
    {
        // No Cuda-aware MPI yet
        auto send_buffer = ddc::create_mirror_view_and_copy(send_span);
        auto recv_buffer = ddc::create_mirror_view(recv_span);
        MPI_Alltoall(
                send_buffer.data_handle(),
                send_buffer.size() / m_comm_size,
                MPI_type_descriptor_t<ElementType>,
                recv_buffer.data_handle(),
                recv_buffer.size() / m_comm_size,
                MPI_type_descriptor_t<ElementType>,
                IMPITranspose<Layout1, Layout2>::m_comm);
        if constexpr (!ddc::is_borrowed_chunk_v<decltype(recv_buffer)>) {
            ddc::parallel_deepcopy(execution_space, recv_span, recv_buffer);
        }
    }

    template <class... DistributedDims>
    ddc::DiscreteDomain<MPIDim<DistributedDims>...> get_distribution(
            ddc::DiscreteDomain<DistributedDims...> local_domain,
            ddc::DiscreteDomain<DistributedDims...> global_domain)
    {
        ddc::DiscreteElement<MPIDim<DistributedDims>...> start(
                ddc::DiscreteElement<MPIDim<DistributedDims>> {0}...);
        ddc::DiscreteVector<MPIDim<DistributedDims>...> size(
                ddc::DiscreteVector<MPIDim<DistributedDims>> {
                        ddc::select<DistributedDims>(global_domain).size()
                        / ddc::select<DistributedDims>(local_domain).size()}...);
        return ddc::DiscreteDomain<MPIDim<DistributedDims>...>(start, size);
    }
};
