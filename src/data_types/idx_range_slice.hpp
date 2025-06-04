// SPDX-License-Identifier: MIT
#pragma once

#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"


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
