

# File mpitransposealltoall.hpp

[**File List**](files.md) **>** [**mpi\_parallelisation**](dir_a35b8fd75f8fad0c2619b083ab571e51.md) **>** [**mpitransposealltoall.hpp**](mpitransposealltoall_8hpp.md)

[Go to the documentation of this file](mpitransposealltoall_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once
#include <numeric>

#include <ddc/ddc.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "impitranspose.hpp"
#include "mpilayout.hpp"
#include "mpitools.hpp"
#include "transpose.hpp"

template <class Layout1, class Layout2>
class MPITransposeAllToAll : public IMPITranspose<Layout1, Layout2>
{
public:
    using idx_range_type1 = typename Layout1::discrete_domain_type;
    using idx_range_type2 = typename Layout2::discrete_domain_type;
    using distributed_idx_range_type1 = typename Layout1::distributed_sub_idx_range;
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

public:
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
        distributed_idx_range_type1 distrib_idx_range(global_idx_range_layout_1);
        if (m_comm_size > distrib_idx_range.size()) {
            throw std::runtime_error("The number of MPI ranks is greater than the number that "
                                     "would be used when maximumly distributing the data");
        }
        // Ensure that the load balancing is good
        if (distrib_idx_range.size() % m_comm_size != 0) {
            throw std::runtime_error("The provided index range cannot be split equally over "
                                     "the specified number of MPI ranks.");
        }
        m_local_idx_range_1
                = m_layout_1.distribute_idx_range(global_idx_range_layout_1, m_comm_size, rank);
        m_local_idx_range_2
                = m_layout_2.distribute_idx_range(global_idx_range_layout_2, m_comm_size, rank);
        m_layout_1_mpi_idx_range = get_distribution(
                distributed_idx_range_type1(m_local_idx_range_1),
                distributed_idx_range_type1(global_idx_range_layout_1));
        m_layout_2_mpi_idx_range = get_distribution(
                distributed_idx_range_type2(m_local_idx_range_2),
                distributed_idx_range_type2(global_idx_range_layout_2));
    }

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
    }

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
        gather_mpi_idx_range_type gather_mpi_idx_range;
        scatter_mpi_idx_range_type scatter_mpi_idx_range;
        if constexpr (std::is_same_v<InLayout, Layout1>) {
            gather_mpi_idx_range = m_layout_1_mpi_idx_range;
            scatter_mpi_idx_range = m_layout_2_mpi_idx_range;
        } else {
            gather_mpi_idx_range = m_layout_2_mpi_idx_range;
            scatter_mpi_idx_range = m_layout_1_mpi_idx_range;
        }

        // Collect the useful subindex ranges described in the fields
        batch_idx_range_type batch_idx_range(get_idx_range(send_field));
        scatter_idx_range_type scatter_idx_range(get_idx_range(recv_field));
        gather_idx_range_type gather_idx_range(get_idx_range(send_field));

        // Build the index ranges describing the function inputs but including the MPI rank information
        input_mpi_idx_range_type input_mpi_idx_range(
                scatter_mpi_idx_range,
                scatter_idx_range,
                gather_idx_range,
                batch_idx_range);
        output_mpi_idx_range_type output_mpi_idx_range(
                gather_mpi_idx_range,
                scatter_idx_range,
                gather_idx_range,
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
                get_const_field(alltoall_send_buffer));

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
    template <
            class ElementType,
            class MPIRecvIdxRange,
            class MPISendIdxRange,
            class MemSpace,
            class ExecSpace>
    void call_all_to_all(
            ExecSpace const& execution_space,
            Field<ElementType, MPIRecvIdxRange, MemSpace> recv_field,
            ConstField<ElementType, MPISendIdxRange, MemSpace> send_field) const
    {
        // No Cuda-aware MPI yet
        auto send_buffer = ddc::create_mirror_view_and_copy(send_field);
        auto recv_buffer = ddc::create_mirror_view(recv_field);
        MPI_Alltoall(
                send_buffer.data_handle(),
                send_buffer.size() / m_comm_size,
                MPI_type_descriptor_t<ElementType>,
                recv_buffer.data_handle(),
                recv_buffer.size() / m_comm_size,
                MPI_type_descriptor_t<ElementType>,
                IMPITranspose<Layout1, Layout2>::m_comm);
        if constexpr (!ddc::is_borrowed_chunk_v<decltype(recv_buffer)>) {
            ddc::parallel_deepcopy(execution_space, recv_field, recv_buffer);
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
```


