// SPDX-License-Identifier: MIT
#pragma once

#include <ddc/kernels/splines.hpp>

#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "derivative_field.hpp"
#include "derivative_field_mem.hpp"
#include "vector_field.hpp"
#include "vector_field_mem.hpp"
#include "vector_index_tools.hpp"

/*
    List of aliases useful to store data in MultipatchType or MultipatchFieldOf.
    Fields and vector fields on each patch should rather be stored using the aliases
    MultipatchField, MultipatchFieldMem, MultipatchVectorField, MultipatchVectorFieldMem, etc.
    with a MultipatchType of the index range aliases below.
*/


// GRIDS -----------------------------------------------------------------------------------------

/// @brief The grid on the first of the Patch's logical dimensions.
template <class Patch>
using Grid1OnPatch = typename Patch::Grid1;

/// @brief The grid on the second of the Patch's logical dimensions.
template <class Patch>
using Grid2OnPatch = typename Patch::Grid2;


// FIELDS ----------------------------------------------------------------------------------------

// --- ON DEVICE ---
/// @brief A Field of doubles defined on the Patch's 2D logical domain.
template <class Patch>
using DFieldOnPatch = DField<typename Patch::IdxRange12>;


// --- ON HOST ---
/// @brief A Field of doubles defined on host on the Patch's 2D logical domain.
template <class Patch>
using DFieldOnPatch_host = host_t<DFieldOnPatch<Patch>>;


// IDX, IDXRANGE ---------------------------------------------------------------------------------

/// @brief An index range over the grids on the Patch's 2D logical domain.
template <class Patch>
using IdxRangeOnPatch = typename Patch::IdxRange12;

/// @brief An index range over the grids on the first Patch's 1D logical domain.
template <class Patch>
using IdxRange1OnPatch = typename Patch::IdxRange1;

/// @brief An index for the grid on the first of the Patch's logical dimensions.
template <class Patch>
using Idx1OnPatch = typename Patch::Idx1;

/// @brief An index range over the B-splines on the Patch's 2D logical domain.
template <class Patch>
using IdxRangeBSOnPatch = typename Patch::IdxRangeBS12;

/**
 * @brief The index range on (dx, x, dy, y) on which a DerivField is defined for the grids on the
 * Patch's 2D logical domain.
 */
template <class Patch>
using IdxRangeDerivOnPatch = IdxRange<
        ddc::Deriv<typename Patch::Dim1>,
        typename Patch::Grid1,
        ddc::Deriv<typename Patch::Dim2>,
        typename Patch::Grid2>;

/**
 * @brief An index range over the derivatives along the Patch's first logical dimension and the
 * grid on the Patch's second logical dimension.
 */
template <class Patch>
using IdxRangeDeriv1Grid2OnPatch
        = IdxRange<ddc::Deriv<typename Patch::Dim1>, typename Patch::Grid2>;

/**
 * @brief An index range over the grid on the Patch's first logical dimension and the derivatives
 * along the Patch's second logical dimension.
 */
template <class Patch>
using IdxRangeGrid1Deriv2OnPatch
        = IdxRange<typename Patch::Grid1, ddc::Deriv<typename Patch::Dim2>>;

/**
 * @brief An index range over the derivatives along both of the Patch's logical dimensions
 * (e.g. for the cross-derivatives).
 */
template <class Patch>
using IdxRangeDeriv12OnPatch
        = IdxRange<ddc::Deriv<typename Patch::Dim1>, ddc::Deriv<typename Patch::Dim2>>;


// IDXRANGESLICE ---------------------------------------------------------------------------------

/// @brief An index range slice over the grids on the first Patch's 1D logical domain.
template <class Patch>
using IdxRange1SliceOnPatch = IdxRangeSlice<typename Patch::Grid1>;

/// @brief An index range slice over the grids on the second Patch's 1D logical domain.
template <class Patch>
using IdxRange2SliceOnPatch = IdxRangeSlice<typename Patch::Grid2>;


// SPLINES ---------------------------------------------------------------------------------------

/// @brief The BSplines dimension on the first of the Patch's logical dimensions.
template <class Patch>
using BSplines1OnPatch = typename Patch::BSplines1;

/// @brief The BSplines dimension on the second of the Patch's logical dimensions.
template <class Patch>
using BSplines2OnPatch = typename Patch::BSplines2;


// FIELDS WITH DERIVATIVES -----------------------------------------------------------------------

/**
 * @brief A field with derivatives containing the function values on the 2D grid, and first derivatives
 * on the boundaries of the 2D grid and the cross-derivatives on the corner of the 2D grid for the 2D grid
 * defined on the given Patch. The data are stored on an index range on (dx, x, dy, y) for a grid on (x,y).
 */
template <class Patch>
using DerivFieldOnPatch = DerivField<double, IdxRangeDerivOnPatch<Patch>>;
