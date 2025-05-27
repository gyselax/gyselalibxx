// SPDX-License-Identifier: MIT
#pragma once

#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"


template <class... DDims>
using IdxRangeSlice = ddc::StridedDiscreteDomain<DDims...>;

/**
 * @brief Check if all elements of the specified IdxRange are found in an index range slice.
 *
 * @param idx_range_slice An index range slice.
 * @param idx_range The index range which may or may not be in this index range slice.
 *
 * @returns bool True if the all elements of the specified IdxRange are found in this index range slice.
 *               False otherwise.
 */
template <class... Dims, class... DDims>
KOKKOS_FUNCTION bool contains(
        IdxRangeSlice<Dims...> const& idx_range_slice,
        IdxRange<DDims...> idx_range)
{
    using tags_seq = ddc::detail::TypeSeq<Dims...>;
    static_assert(
            (ddc::in_tags_v<DDims, tags_seq> && ...),
            "Requested Tag absent from IdxRangeSlice");
    return IdxRangeSlice<DDims...>(idx_range_slice).contains(idx_range.front())
           && (((idx_range.template extent<DDims>().value() == 1)
                || (ddc::select<DDims>(idx_range_slice.strides()) == 1))
               && ...);
}

/// A class to create a IdxRangeSlice type from a TypeSeq.
template <class T>
struct IdxRangeToSlice;

template <class... Dims>
struct IdxRangeToSlice<ddc::detail::TypeSeq<Dims...>>
{
    using value = IdxRangeSlice<Dims...>;
};

/// @brief A templated type for creating a IdxRangeSlice from a TypeSeq.
template <class T>
using to_subidx_range_collection = typename IdxRangeToSlice<T>::value;
