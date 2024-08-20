// SPDX-License-Identifier: MIT

#pragma once
#include "ddc_aliases.hpp"


/// @brief Type for MultipatchType: A Field defined on the Patch's 2D logical domain.
template <class Patch>
using DFieldOnPatch = DField<typename Patch::IdxRange12>;

/// @brief Type for MultipatchType: A constant Field defined on the Patch's 2D logical domain.
template <class Patch>
using DConstFieldOnPatch = DConstField<typename Patch::IdxRange12>;

/// @brief Type for MultipatchType: An index range over the grids on the Patch's 2D logical domain.
template <class Patch>
using IdxRangeOnPatch = typename Patch::IdxRange12;

/// @brief Type for MultipatchType: An index for the grid on the first of the Patch's logical dimensions.
template <class Patch>
using Idx1OnPatch = typename Patch::Idx1;

/// @brief Type for MultipatchType: The BSplines on the first of the Patch's logical dimensions.
template <class Patch>
using BSplines1OnPatch = typename Patch::BSplines1;

/// @brief Type for MultipatchType: The grid on the first of the Patch's logical dimensions.
template <class Patch>
using Grid1OnPatch = typename Patch::Grid1;

/// @brief Type for MultipatchType: A constant Field of doubles on the first of the Patch's logical dimensions.
template <class Patch>
using DConstField1OnPatch = DConstField<typename Patch::IdxRange1>;

/**
 * @brief Type for MultipatchType: A field of spline coefficients batched over the second of the Patch's
 * logical dimensions for a spline defined on the first of the Patch's logical dimensions.
 */
template <class Patch>
using SplineCoeff1OnPatch_2D = DField<IdxRange<typename Patch::BSplines1, typename Patch::Grid2>>;

/**
 * @brief Type for MultipatchType: A field of spline coefficients for a non-batched spline defined
 * on the first of the Patch's logical dimensions.
 */
template <class Patch>
using SplineCoeff1OnPatch_1D = DField<typename Patch::IdxRangeBS1>;
