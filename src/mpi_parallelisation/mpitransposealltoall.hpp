// SPDX-License-Identifier: MIT
#pragma once
#include <numeric>
#include <vector>

#include <ddc/ddc.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "impitranspose.hpp"
#include "mpilayout.hpp"
#include "mpitools.hpp"
#include "transpose.hpp"

/**
 * @brief A class describing an operator for converting from/to different MPI layouts using AlltoAll.
 *
 * This class implements a basic AlltoAll operator and currently only works with a basic MPILayout.
 * There may be more MPI processes than points in the distributed dimensions of one (or both)
 * of the layouts, in which case some processes are inactive (they hold no data). Such processes
 * must still call the transpose operator (it is collective) but they send/receive no data
 * in the corresponding layout.
 *
 * @tparam Layout1 One of the MPI layouts.
 * @tparam Layout2 The other MPI layouts.
 */
template <class Layout1, class Layout2>
class MPITransposeAllToAll : public IMPITranspose<Layout1, Layout2>
{
public:
    /// The type of the index range of the first MPI layout.
    using idx_range_type1 = typename Layout1::discrete_domain_type;
    /// The type of the index range of the second MPI layout.
    using idx_range_type2 = typename Layout2::discrete_domain_type;
    /// The type of the index range of the first MPI layout.
    using distributed_idx_range_type1 = typename Layout1::distributed_sub_idx_range;
    /// The type of the index range of the second MPI layout.
    using distributed_idx_range_type2 = typename Layout2::distributed_sub_idx_range;

private:
    using layout_1_mpi_dims = ddcHelper::
            apply_template_to_type_seq_t<MPIDim, typename Layout1::distributed_type_seq>;
    using layout_2_mpi_dims = ddcHelper::
            apply_template_to_type_seq_t<MPIDim, typename Layout2::distributed_type_seq>;
    using layout_1_mpi_idx_range_type
            = ddc::detail::convert_type_seq_to_discrete_domain_t<layout_1_mpi_dims>;
    using layout_2_mpi_idx_range_type
            = ddc::detail::convert_type_seq_to_discrete_domain_t<layout_2_mpi_dims>;

private:
    int m_comm_size;
    Layout1 m_layout_1;
    Layout2 m_layout_2;
    idx_range_type1 m_local_idx_range_1;
    idx_range_type2 m_local_idx_range_2;
    layout_1_mpi_idx_range_type m_layout_1_mpi_idx_range;
    layout_2_mpi_idx_range_type m_layout_2_mpi_idx_range;
    // The distributed index range handled by an active process (rank 0 used as reference).
    distributed_idx_range_type1 m_nonempty_distributed_idx_range_1;
    distributed_idx_range_type2 m_nonempty_distributed_idx_range_2;

public:
    /**
     * @brief A constructor for the transpose operator.
     *
     * @param global_idx_range The global index range of the data across processes provided in either layout.
     * @param comm The MPI communicator.
     */
    template <class IdxRange>
    MPITransposeAllToAll(IdxRange global_idx_range, MPI_Comm comm)
        : IMPITranspose<Layout1, Layout2>(comm)
    {
        static_assert(
                std::is_same_v<
                        IdxRange,
                        idx_range_type1> || std::is_same_v<IdxRange, idx_range_type2>,
                "The initialisation global idx_range should be described by one of the layouts");
        idx_range_type1 global_idx_range_layout_1(global_idx_range);
        idx_range_type2 global_idx_range_layout_2(global_idx_range);
        int rank;
        MPI_Comm_size(comm, &m_comm_size);
        MPI_Comm_rank(comm, &rank);
        // Index range of held by this MPI rank (may be empty)
        m_local_idx_range_1
                = m_layout_1.distribute_idx_range(global_idx_range_layout_1, m_comm_size, rank);
        m_local_idx_range_2
                = m_layout_2.distribute_idx_range(global_idx_range_layout_2, m_comm_size, rank);
        // Reference, non-empty index ranges, based on rank 0
        m_nonempty_distributed_idx_range_1 = distributed_idx_range_type1(
                m_layout_1.distribute_idx_range(global_idx_range_layout_1, m_comm_size, 0));
        m_nonempty_distributed_idx_range_2 = distributed_idx_range_type2(
                m_layout_2.distribute_idx_range(global_idx_range_layout_2, m_comm_size, 0));
        // Distributed index ranges
        m_layout_1_mpi_idx_range = get_distribution(
                m_nonempty_distributed_idx_range_1,
                distributed_idx_range_type1(global_idx_range_layout_1));
        m_layout_2_mpi_idx_range = get_distribution(
                m_nonempty_distributed_idx_range_2,
                distributed_idx_range_type2(global_idx_range_layout_2));

        assert(m_layout_1_mpi_idx_range.size() == m_layout_1.n_active_ranks());
        assert(m_layout_2_mpi_idx_range.size() == m_layout_2.n_active_ranks());
    }

    /**
     * @brief Getter for the local index range.
     *
     * @tparam Layout The layout whose index range should be retrieved.
     *
     * @returns The local index range in the specified MPI layout.
     */
    template <class Layout>
    auto get_local_idx_range() const
    {
        static_assert(
                std::is_same_v<Layout, Layout1> || std::is_same_v<Layout, Layout2>,
                "Transpose class does not handle requested layout");
        if constexpr (std::is_same_v<Layout, Layout1>) {
            return m_local_idx_range_1;
        } else if constexpr (std::is_same_v<Layout, Layout2>) {
            return m_local_idx_range_2;
        }
    }

    /**
     * @brief An operator which transposes from one layout to another.
     *
     * If the dimensions are ordered differently in the index ranges of the two layouts
     * then this function can be used. If the index ranges have the same type then the
     * transpose_to function must be used to disambiguate.
     *
     * @param[in] execution_space The execution space (Host/Device) where the code will run.
     * @param[out] recv_field The chunk which will describe the data in the new layout. This
     *                      data is gathered from the MPI processes.
     * @param[in] send_field The chunk describing the data in the current layout. This data
     *                      will be scattered to other MPI processes.
     */
    template <
            class ElementType,
            class InIdxRange,
            class IdxRangeOut,
            class MemSpace,
            class ExecSpace>
    void operator()(
            ExecSpace const& execution_space,
            Field<ElementType, IdxRangeOut, MemSpace> recv_field,
            ConstField<ElementType, InIdxRange, MemSpace> send_field) const
    {
        Kokkos::Profiling::pushRegion("(GSLX) MpiTranspose");
        static_assert(!std::is_same_v<InIdxRange, IdxRangeOut>);
        static_assert(
                (std::is_same_v<InIdxRange, typename Layout1::discrete_domain_type>)
                || (std::is_same_v<InIdxRange, typename Layout2::discrete_domain_type>));
        static_assert(
                (std::is_same_v<IdxRangeOut, typename Layout1::discrete_domain_type>)
                || (std::is_same_v<IdxRangeOut, typename Layout2::discrete_domain_type>));
        using OutLayout = std::conditional_t<
                std::is_same_v<IdxRangeOut, typename Layout1::discrete_domain_type>,
                Layout1,
                Layout2>;
        this->transpose_to<OutLayout>(execution_space, recv_field, send_field);
        Kokkos::Profiling::popRegion();
    }

    /**
     * @brief An operator which transposes from one layout to another.
     *
     * If the dimensions are ordered differently in the index ranges of the two layouts
     * then this function can be used. If the index ranges have the same type then the
     * transpose_to function must be used to disambiguate.
     *
     * @tparam OutLayout The layout that the data should be transposed to.
     *
     * @param[in] execution_space The execution space (Host/Device) where the code will run.
     * @param[out] recv_field The chunk which will describe the data in the new layout. This
     *                      data is gathered from the MPI processes.
     * @param[in] send_field The chunk describing the data in the current layout. This data
     *                      will be scattered to other MPI processes.
     */
    template <class OutLayout, class ElementType, class MemSpace, class ExecSpace, class InIdxRange>
    void transpose_to(
            ExecSpace const& execution_space,
            Field<ElementType, typename OutLayout::discrete_domain_type, MemSpace> recv_field,
            ConstField<ElementType, InIdxRange, MemSpace> send_field) const
    {
        using InLayout = std::conditional_t<std::is_same_v<OutLayout, Layout1>, Layout2, Layout1>;
        /*****************************************************************
         * Validate input
         *****************************************************************/
        static_assert((std::is_same_v<InLayout, Layout1>) || (std::is_same_v<InLayout, Layout2>));
        static_assert((std::is_same_v<OutLayout, Layout1>) || (std::is_same_v<OutLayout, Layout2>));
        static_assert(Kokkos::SpaceAccessibility<ExecSpace, MemSpace>::accessible);

        static_assert(std::is_same_v<InIdxRange, typename InLayout::discrete_domain_type>);
        using IdxRangeOut = typename OutLayout::discrete_domain_type;

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
        // Get complete set of tags from index ranges
        using input_ordered_dims = ddc::to_type_seq_t<InIdxRange>;
        using output_ordered_dims = ddc::to_type_seq_t<IdxRangeOut>;
        // Find tags that are neither scattered nor gathered
        using batch_dims = ddc::type_seq_remove_t<
                ddc::type_seq_remove_t<input_ordered_dims, dims_to_scatter>,
                dims_to_gather>;
        // Insert the MPI tags into the received data to introduce the dimension split
        using input_mpi_idx_range_tags
                = insert_mpi_tags_into_seq_t<scatter_mpi_dims, input_ordered_dims>;
        using output_mpi_idx_range_tags
                = insert_mpi_tags_into_seq_t<gather_mpi_dims, output_ordered_dims>;
        // During the alltoall call the MPI tags must be first to correctly send/receive data
        // but data does not reorder itself
        using input_alltoall_dim_order
                = ddc::type_seq_merge_t<scatter_mpi_dims, input_ordered_dims>;
        using output_alltoall_dim_order
                = ddc::type_seq_merge_t<gather_mpi_dims, input_ordered_dims>;

        /*****************************************************************
         * Define useful index range types
         *****************************************************************/
        // Get the index ranges of objects which are distributed across processes
        using gather_idx_range_type = typename InLayout::distributed_sub_idx_range;
        // Get the index ranges of objects which will be distributed across processes
        using scatter_idx_range_type = typename OutLayout::distributed_sub_idx_range;
        // Get the index range of objects that are not distributed across processes
        using batch_idx_range_type = ddc::detail::convert_type_seq_to_discrete_domain_t<batch_dims>;
        using gather_mpi_idx_range_type
                = ddc::detail::convert_type_seq_to_discrete_domain_t<gather_mpi_dims>;
        using scatter_mpi_idx_range_type
                = ddc::detail::convert_type_seq_to_discrete_domain_t<scatter_mpi_dims>;
        // Get the index range containing MPI tags which can be used to describe the function input
        using input_mpi_idx_range_type
                = ddc::detail::convert_type_seq_to_discrete_domain_t<input_mpi_idx_range_tags>;
        // Get the index range containing MPI tags which can be used to describe the function output
        using output_mpi_idx_range_type
                = ddc::detail::convert_type_seq_to_discrete_domain_t<output_mpi_idx_range_tags>;
        // Get the index range that should be used when scattering data
        using input_alltoall_idx_range_type
                = ddc::detail::convert_type_seq_to_discrete_domain_t<input_alltoall_dim_order>;
        // Get the index range that should be used when gathering data
        using output_alltoall_idx_range_type
                = ddc::detail::convert_type_seq_to_discrete_domain_t<output_alltoall_dim_order>;

        /*****************************************************************
         * Build index ranges
         *****************************************************************/
        // Collect the artificial dimension describing the MPI rank where the scattered information
        // will be sent to or where the gathered information will be collected from
        // Only the active ranks appear in these index ranges.
        gather_mpi_idx_range_type gather_mpi_idx_range;
        scatter_mpi_idx_range_type scatter_mpi_idx_range;
        // Index ranges handled by an active rank (non-empty even if this rank is inactive).
        gather_idx_range_type nonempty_gather_idx_range;
        scatter_idx_range_type nonempty_scatter_idx_range;

        if constexpr (std::is_same_v<InLayout, Layout1>) {
            gather_mpi_idx_range = m_layout_1_mpi_idx_range;
            scatter_mpi_idx_range = m_layout_2_mpi_idx_range;
            nonempty_gather_idx_range = m_nonempty_distributed_idx_range_1;
            nonempty_scatter_idx_range = m_nonempty_distributed_idx_range_2;
        } else {
            gather_mpi_idx_range = m_layout_2_mpi_idx_range;
            scatter_mpi_idx_range = m_layout_1_mpi_idx_range;
            nonempty_gather_idx_range = m_nonempty_distributed_idx_range_2;
            nonempty_scatter_idx_range = m_nonempty_distributed_idx_range_1;
        }

        // Collect the useful subindex ranges described in the fields.
        batch_idx_range_type batch_idx_range(get_idx_range(send_field));
        scatter_idx_range_type local_scatter_idx_range(get_idx_range(recv_field));
        gather_idx_range_type local_gather_idx_range(get_idx_range(send_field));

        // Build the index ranges describing the function inputs but including the MPI rank information.
        // The data sent to each active rank in the output layout is described by the reference,
        // non-empty scatter index range.
        input_mpi_idx_range_type input_mpi_idx_range(
                scatter_mpi_idx_range,
                nonempty_scatter_idx_range,
                local_gather_idx_range,
                batch_idx_range);
        output_mpi_idx_range_type output_mpi_idx_range(
                gather_mpi_idx_range,
                local_scatter_idx_range,
                nonempty_gather_idx_range,
                batch_idx_range);
        assert(input_mpi_idx_range.size() == send_field.size());
        assert(output_mpi_idx_range.size() == recv_field.size());

        // Create the index ranges used during the alltoall call (the MPI rank index range is first in this layout)
        input_alltoall_idx_range_type input_alltoall_idx_range(input_mpi_idx_range);
        output_alltoall_idx_range_type output_alltoall_idx_range(output_mpi_idx_range);

        /*****************************************************************
         * Create views on the function inputs with the index ranges including the MPI rank information
         *****************************************************************/
        ConstField<ElementType, input_mpi_idx_range_type, MemSpace>
                send_mpi_field(send_field.data_handle(), input_mpi_idx_range);
        Field<ElementType, output_mpi_idx_range_type, MemSpace>
                recv_mpi_field(recv_field.data_handle(), output_mpi_idx_range);

        /*****************************************************************
         * Transpose data (both on the rank and between ranks)
         *****************************************************************/
        // Create views or copies of the function inputs such that they are laid out on the
        // index range used during the alltoall call
        auto alltoall_send_buffer = ddcHelper::create_transpose_mirror_view_and_copy<
                input_alltoall_idx_range_type>(execution_space, send_mpi_field);
        auto alltoall_recv_buffer = ddcHelper::create_transpose_mirror<
                output_alltoall_idx_range_type>(execution_space, recv_mpi_field);

        // Call the MPI AlltoAll routine
        call_all_to_all(
                execution_space,
                get_field(alltoall_recv_buffer),
                get_const_field(alltoall_send_buffer),
                scatter_mpi_idx_range.size(),
                gather_mpi_idx_range.size());

        // If alltoall_recv_buffer owns its data (not a view) then copy the results back to
        // recv_mpi_field which is a view on recv_field, the function output
        if constexpr (!ddc::is_borrowed_chunk_v<decltype(alltoall_recv_buffer)>) {
            transpose_layout(
                    execution_space,
                    recv_mpi_field,
                    get_const_field(alltoall_recv_buffer));
        }
    }

private:
    /**
     * @brief Function handling the MPI call.
     *
     * @param[in] execution_space The execution space (Host/Device) where the code will run.
     * @param[out] recv_field The chunk which will describe the data in the new layout. This
     *                      data is gathered from the MPI processes.
     * @param[in] send_field The chunk describing the data in the current layout. This data
     *                      will be scattered to other MPI processes.
     * @param[in] n_send_ranks Number of active ranks sending data.
     * @param[in] n_recv_ranks Number of active ranks receiving data.
     */
    template <
            class ElementType,
            class MPIRecvIdxRange,
            class MPISendIdxRange,
            class MemSpace,
            class ExecSpace>
    void call_all_to_all(
            ExecSpace const& execution_space,
            Field<ElementType, MPIRecvIdxRange, MemSpace> recv_field,
            ConstField<ElementType, MPISendIdxRange, MemSpace> send_field,
            int n_send_ranks,
            int n_recv_ranks) const
    {
        Kokkos::fence("fencing before mpi_all_to_all");

        if (n_send_ranks == m_comm_size && n_recv_ranks == m_comm_size) {
            // All ranks involved are active so MPI_Alltoall can be used
            MPI_Alltoall(
                    send_field.data_handle(),
                    send_field.size() / m_comm_size,
                    MPI_type_descriptor_t<ElementType>,
                    recv_field.data_handle(),
                    recv_field.size() / m_comm_size,
                    MPI_type_descriptor_t<ElementType>,
                    IMPITranspose<Layout1, Layout2>::m_comm);
        } else {
            // Some ranks are inactive, so MPI_Alltoallv must be used
            int const send_block = send_field.size() / n_send_ranks;
            int const recv_block = recv_field.size() / n_recv_ranks;
            std::vector<int> send_counts(m_comm_size, 0);
            std::vector<int> send_displs(m_comm_size, 0);
            std::vector<int> recv_counts(m_comm_size, 0);
            std::vector<int> recv_displs(m_comm_size, 0);
            for (int i = 0; i < n_send_ranks; ++i) {
                send_counts[i] = send_block;
                send_displs[i] = i * send_block;
            }
            for (int i = 0; i < n_recv_ranks; ++i) {
                recv_counts[i] = recv_block;
                recv_displs[i] = i * recv_block;
            }
            MPI_Alltoallv(
                    send_field.data_handle(),
                    send_counts.data(),
                    send_displs.data(),
                    MPI_type_descriptor_t<ElementType>,
                    recv_field.data_handle(),
                    recv_counts.data(),
                    recv_displs.data(),
                    MPI_type_descriptor_t<ElementType>,
                    IMPITranspose<Layout1, Layout2>::m_comm);
        }
    }

    template <class... DistributedDims>
    IdxRange<MPIDim<DistributedDims>...> get_distribution(
            IdxRange<DistributedDims...> local_idx_range,
            IdxRange<DistributedDims...> global_idx_range) const
    {
        Idx<MPIDim<DistributedDims>...> start(Idx<MPIDim<DistributedDims>> {0}...);
        IdxStep<MPIDim<DistributedDims>...> size(IdxStep<MPIDim<DistributedDims>> {
                ddc::select<DistributedDims>(global_idx_range).size()
                / ddc::select<DistributedDims>(local_idx_range).size()}...);
        return IdxRange<MPIDim<DistributedDims>...>(start, size);
    }
};
