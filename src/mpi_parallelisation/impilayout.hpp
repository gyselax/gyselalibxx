// SPDX-License-Identifier: MIT
#pragma once
#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"

/**
 * @brief A super class describing a way in which data may be laid out across MPI processes.
 *
 * @tparam DataIdxRange The index range on which the data is defined.
 * @tparam DistributedDim The tags of the discrete dimensions which are distributed
 *              across MPI processes.
 */
template <class DataIdxRange, class... DistributedDim>
class IMPILayout
{
    static_assert(ddc::is_discrete_domain_v<DataIdxRange>);

public:
    /// The index range of the data
    using discrete_domain_type = DataIdxRange;
    /// The index range of the distributed section of the data
    using distributed_sub_idx_range = IdxRange<DistributedDim...>;
    /// A type sequence describing the dimensions which are distributed across MPI processes.
    using distributed_type_seq = ddc::detail::TypeSeq<DistributedDim...>;

private:
    using type_seq = typename ddc::to_type_seq_t<discrete_domain_type>;

    template <std::size_t... I>
    static constexpr bool check_distributed_idx_ranges_are_first(
            std::integer_sequence<std::size_t, I...>)
    {
        return (ddc::in_tags_v<ddc::type_seq_element_t<I, type_seq>, distributed_type_seq> && ...);
    }

public:
    /// The number of dimensions that are distributed across MPI processes.
    static constexpr int n_distributed_dimensions = sizeof...(DistributedDim);
    /**
     * @brief A flag to indicate whether the distributed dimensions are the dimensions which
     * are the furthest from being contiguous in memory.
     */
    static constexpr bool distributed_idx_ranges_are_first = check_distributed_idx_ranges_are_first(
            std::make_index_sequence<n_distributed_dimensions>());
};
