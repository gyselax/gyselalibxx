// SPDX-License-Identifier: MIT

#pragma once
#include "ddc_aliases.hpp"


/// @brief Type for MultipatchType: tuple of Fields defined on the Patch logical domain.
template <class Patch>
using DFieldOnPatch = DField<typename Patch::IdxRange12>;

/// @brief Type for MultipatchType: tuple of Fields defined on the Patch logical domain.
template <class Patch>
using DConstFieldOnPatch = DConstField<typename Patch::IdxRange12>;


/// @brief Type for MultipatchType: tuple of index ranges defined on the Patch logical domain.
template <class Patch>
using IdxRangeOnPatch = typename Patch::IdxRange12;

/// @brief Type for MultipatchType: tuple of indexes on Dim1 defined on the Patch logical domain.
template <class Patch>
using Idx1OnPatch = typename Patch::Idx1;
