

# File mpilayout.hpp

[**File List**](files.md) **>** [**mpi\_parallelisation**](dir_a35b8fd75f8fad0c2619b083ab571e51.md) **>** [**mpilayout.hpp**](mpilayout_8hpp.md)

[Go to the documentation of this file](mpilayout_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

#include <numeric>
#include <sstream>

#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"
#include "impilayout.hpp"

template <class IdxRangeData, class... DistributedDim>
class MPILayout : public IMPILayout<IdxRangeData, DistributedDim...>
{
    static_assert(ddc::is_discrete_domain_v<IdxRangeData>);

    using base_type = IMPILayout<IdxRangeData, DistributedDim...>;

public:
    using idx_range_type = IdxRangeData;
    using distributed_sub_idx_range = typename base_type::distributed_sub_idx_range;
    using distributed_type_seq = typename base_type::distributed_type_seq;

public:
    static idx_range_type distribute_idx_range(
            idx_range_type global_idx_range,
            int comm_size,
            int rank)
    {
        distributed_sub_idx_range distrib_idx_range(global_idx_range);
        if (distrib_idx_range.size() % comm_size != 0) {
            std::ostringstream error_msg;
            error_msg << "The provided index range cannot be split equally over the specified "
                         "number of MPI ranks ("
                      << distrib_idx_range.extents() << " is not divisible by " << comm_size << ")";
            throw std::runtime_error(error_msg.str());
        }
        return internal_distribute_idx_range(global_idx_range, comm_size, rank);
    }



protected:
    template <class HeadTag>
    static IdxRange<HeadTag> internal_distribute_idx_range(
            IdxRange<HeadTag> global_idx_range,
            int comm_size,
            int rank)
    {
        if constexpr (ddc::in_tags_v<HeadTag, distributed_type_seq>) {
            if (global_idx_range.size() % comm_size != 0) {
                throw std::runtime_error("The provided index range cannot be split equally over "
                                         "the specified number of MPI ranks.");
            }
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

    template <class HeadTag, class... Tags, std::enable_if_t<(sizeof...(Tags) > 0), bool> = true>
    static IdxRange<HeadTag, Tags...> internal_distribute_idx_range(
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
```


