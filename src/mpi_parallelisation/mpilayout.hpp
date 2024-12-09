// SPDX-License-Identifier: MIT
#pragma once

#include <numeric>

#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"
#include "impilayout.hpp"

/**
 * @brief A class describing a way in which data may be laid out across MPI processes.
 *
 * This class describes the simplest way of laying out data. In this layout the data is
 * distributed along dimensions in order until the data is distributed. It is possible that
 * the data may not be distributed along some of the requested dimensions if such distribution
 * is not necessary.
 * For example if we distribute dimensions R and Theta of a (R, Theta, Phi, Vpar, Mu) grid
 * of size 256 x 1024 x 256 x 64 x 8 over 64 processes, the theta dimension will not be
 * distributed.
 *
 * The data is distributed such that it is maximally distributed over each dimension
 * (in order) such that the size of each local index range along that dimension is the same for
 * all processes.
 * For example if we distribute dimensions X and Y of a (X, Y, Z) grid of size (10, 15, 4)
 * over 6 processes, the X dimension will be distributed over 2 processes and the Y
 * dimension will be distributed over 3 processes.
 *
 * @tparam IdxRangeData The IdxRange on which the data is defined.
 * @tparam DistributedDim The tags of the discrete dimensions which are distributed
 *              across MPI processes.
 */
template <class IdxRangeData, class... DistributedDim>
class MPILayout : public IMPILayout<IdxRangeData, DistributedDim...>
{
    static_assert(ddc::is_discrete_domain_v<IdxRangeData>);

    using base_type = IMPILayout<IdxRangeData, DistributedDim...>;

public:
    /// The index range of the data
    using idx_range_type = IdxRangeData;
    /// The index range of the distributed section of the data
    using distributed_sub_idx_range = typename base_type::distributed_sub_idx_range;
    /// A type sequence describing the dimensions which are distributed across MPI processes.
    using distributed_type_seq = typename base_type::distributed_type_seq;

public:
    /**
     * @brief Get the distributed index range which follows the chosen layout.
     *
     * @param[in] global_idx_range The global (non-distributed) index range.
     * @param[in] comm_size The number of MPI processes over which the data is distributed.
     * @param[in] rank The rank of the current MPI process.
     *
     * @returns The distributed index range.
     */
    idx_range_type distribute_idx_range(idx_range_type global_idx_range, int comm_size, int rank)
    {
        return internal_distribute_idx_range(global_idx_range, comm_size, rank);
    }



protected:
    /**
     * @brief Distribute a 1D index range over the MPI processes.
     *
     * @param[in] global_idx_range The index range to be distributed.
     * @param[in] comm_size The number of processes over which the data should be disctributed.
     * @param[in] rank The rank of the process within the current group of processes
     *
     * @returns The distributed index range.
     */
    template <class HeadTag>
    IdxRange<HeadTag> internal_distribute_idx_range(
            IdxRange<HeadTag> global_idx_range,
            int comm_size,
            int rank)
    {
        if constexpr (ddc::in_tags_v<HeadTag, distributed_type_seq>) {
            assert(global_idx_range.size() % comm_size == 0);
            IdxStep<HeadTag> elems_on_dim(global_idx_range.size() / comm_size);
            IdxRange<HeadTag>
                    local_idx_range(global_idx_range.front() + rank * elems_on_dim, elems_on_dim);
            return local_idx_range;
        } else {
            // Data is not actually distributed as it handles the case of an index range which is not defined on a distributed dimension.
            assert(comm_size == 1);
            return global_idx_range;
        }
    }

    /**
     * @brief Distribute the index range over the MPI processes.
     *
     * This function is called recursively. At each pass it distributes the index range over the
     * first dimension in the discrete index range. The remaining dimensions and processes are
     * then handled in the recursive call.
     *
     * @param[in] idx_range The index range to be distributed.
     * @param[in] comm_size The number of processes over which the data should be disctributed.
     * @param[in] rank The rank of the process within the current group of processes
     *
     * @returns The distributed index range.
     */
    template <class HeadTag, class... Tags, std::enable_if_t<(sizeof...(Tags) > 0), bool> = true>
    IdxRange<HeadTag, Tags...> internal_distribute_idx_range(
            IdxRange<HeadTag, Tags...> idx_range,
            int comm_size,
            int rank)
    {
        IdxRange<HeadTag> global_idx_range_along_dim = ddc::select<HeadTag>(idx_range);
        IdxRange<HeadTag> local_idx_range_along_dim;
        IdxRange<Tags...> remaining_idx_range;

        if constexpr (ddc::in_tags_v<HeadTag, distributed_type_seq>) {
            // The number of MPI processes along this dimension
            int n_ranks_along_dim = std::gcd(comm_size, global_idx_range_along_dim.size());
            // The number of MPI processes along all subsequent dimensions
            int n_elems_lower_dims = comm_size / n_ranks_along_dim;
            // The rank index for the MPI process along this dimension
            int rank_along_dim = rank / n_elems_lower_dims;
            // The rank index for the MPI process along all subsequent dimensions
            int remaining_rank = rank % n_elems_lower_dims;
            // Calculate the local index range
            IdxStep<HeadTag> elems_on_dim(global_idx_range_along_dim.size() / n_ranks_along_dim);
            Idx<HeadTag> distrib_start(
                    global_idx_range_along_dim.front() + rank_along_dim * elems_on_dim);
            local_idx_range_along_dim = IdxRange<HeadTag>(distrib_start, elems_on_dim);
            // Calculate the index range for the subsequent dimensions
            IdxRange<Tags...> remaining_dims = ddc::select<Tags...>(idx_range);
            remaining_idx_range = internal_distribute_idx_range(
                    remaining_dims,
                    n_elems_lower_dims,
                    remaining_rank);
        } else {
            // Calculate the local index range
            local_idx_range_along_dim = global_idx_range_along_dim;
            // Calculate the index range for the subsequent dimensions
            remaining_idx_range = internal_distribute_idx_range(
                    ddc::select<Tags...>(idx_range),
                    comm_size,
                    rank);
        }
        return IdxRange<HeadTag, Tags...>(local_idx_range_along_dim, remaining_idx_range);
    }
};
