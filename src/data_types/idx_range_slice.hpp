// SPDX-License-Identifier: MIT
#pragma once
#include <map>
#include <vector>

#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"


template <class... DDims>
using IdxRangeSlice = ddc::StridedDiscreteDomain<DDims...>;

template <class QueryDim, class... Dims>
KOKKOS_FUNCTION std::size_t get_index(IdxRangeSlice<Dims...> const& rng, Idx<QueryDim> elem)
{
    assert(rng.contains(elem));
    return IdxRangeSlice<QueryDim>(rng).distance_from_front(elem).value();
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
