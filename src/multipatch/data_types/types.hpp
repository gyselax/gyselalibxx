// SPDX-License-Identifier: MIT

#pragma once
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"


// GRIDS -----------------------------------------------------------------------------------------

/// @brief Type for MultipatchType: The grid on the first of the Patch's logical dimensions.
template <class Patch>
using Grid1OnPatch = typename Patch::Grid1;

/// @brief Type for MultipatchType: The grid on the second of the Patch's logical dimensions.
template <class Patch>
using Grid2OnPatch = typename Patch::Grid2;



// FIELDS ----------------------------------------------------------------------------------------

/// @brief A FieldMem defined on the Patch's 2D logical domain. This type should not be used in a MultipatchType.
template <class Patch>
using DFieldMemOnPatch = DFieldMem<typename Patch::IdxRange12>;

/// @brief Type for MultipatchType: A Field defined on the Patch's 2D logical domain.
template <class Patch>
using DFieldOnPatch = DField<typename Patch::IdxRange12>;

/// @brief Type for MultipatchType: A constant Field defined on the Patch's 2D logical domain.
template <class Patch>
using DConstFieldOnPatch = DConstField<typename Patch::IdxRange12>;

/// @brief Type for MultipatchType: A Field defined on host on the Patch's 2D logical domain.
template <class Patch>
using DFieldOnPatch_host = host_t<DFieldOnPatch<Patch>>;


/// @brief Type for MultipatchType: A Field of doubles on the first of the Patch's logical dimensions.
template <class Patch>
using DField1OnPatch = DField<typename Patch::IdxRange1>;

/// @brief Type for MultipatchType: A constant Field of doubles on the first of the Patch's logical dimensions.
template <class Patch>
using DConstField1OnPatch = DConstField<typename Patch::IdxRange1>;



// IDX, IDXRANGE ---------------------------------------------------------------------------------

/// @brief Type for MultipatchType: An index range over the grids on the Patch's 2D logical domain.
template <class Patch>
using IdxRangeOnPatch = typename Patch::IdxRange12;


/// @brief Type for MultipatchType: An index range over the grids on the first Patch's 1D logical domain.
template <class Patch>
using IdxRange1OnPatch = typename Patch::IdxRange1;

/// @brief Type for MultipatchType: An index for the grid on the first of the Patch's logical dimensions.
template <class Patch>
using Idx1OnPatch = typename Patch::Idx1;



// COORDINATES -----------------------------------------------------------------------------------

/// @brief Type for MultipatchType: Field of 2D coordinates defined on the 2D Patch's logical domain.
template <class Patch>
using CoordFieldOnPatch = Field<typename Patch::Coord12, typename Patch::IdxRange12>;

/// @brief Type for MultipatchType: ConstField of 2D coordinates defined on the 2D Patch's logical domain.
template <class Patch>
using CoordConstFieldOnPatch = ConstField<typename Patch::Coord12, typename Patch::IdxRange12>;


/// @brief Type for MultipatchType: Field of 1D coordinates defined on the defined
/// on domain of the Patch's first logical dimension.
template <class Patch>
using Coord1Field1OnPatch_1D = Field<typename Patch::Coord1, typename Patch::IdxRange1>;



// SPLINES ---------------------------------------------------------------------------------------

/// @brief Type for MultipatchType: The BSplines on the first of the Patch's logical dimensions.
template <class Patch>
using BSplines1OnPatch = typename Patch::BSplines1;

/// @brief Type for MultipatchType: The BSplines on the second of the Patch's logical dimensions.
template <class Patch>
using BSplines2OnPatch = typename Patch::BSplines2;


/**
 * @brief Type for MultipatchType: A field of 2D spline coefficients for a non-batched spline defined
 * on both of the Patch's logical dimensions.
 */
template <class Patch>
using SplineCoeffOnPatch_2D = DField<typename Patch::IdxRangeBS12>;

/**
 * @brief Type for MultipatchType: A field of 2D spline coefficients for a non-batched spline defined
 * on both of the Patch's logical dimensions.
 */
template <class Patch>
using ConstSplineCoeffOnPatch_2D = DConstField<typename Patch::IdxRangeBS12>;

/**
 * @brief Type for MultipatchType: A field of 2D spline coefficients for a non-batched spline defined
 * on both of the Patch's logical dimensions. Defined on host. 
 */
template <class Patch>
using ConstSplineCoeffOnPatch_2D_host = host_t<ConstSplineCoeffOnPatch_2D<Patch>>;

/**
 * @brief Type for MultipatchType: A field of spline coefficients batched over the second of the Patch's
 * logical dimensions for a spline defined on the first of the Patch's logical dimensions.
 */
template <class Patch>
using SplineCoeff1OnPatch_2D = DField<IdxRange<typename Patch::BSplines1, typename Patch::Grid2>>;

/**
 * @brief Type for MultipatchType: A field of spline coefficients batched over the second of the Patch's
 * logical dimensions for a spline defined on the first of the Patch's logical dimensions.
 */
template <class Patch>
using ConstSplineCoeff1OnPatch_2D
        = DConstField<IdxRange<typename Patch::BSplines1, typename Patch::Grid2>>;


/**
 * @brief Type for MultipatchType: A field of spline coefficients for a non-batched spline defined
 * on the first of the Patch's logical dimensions.
 */
template <class Patch>
using SplineCoeff1OnPatch_1D = DField<typename Patch::IdxRangeBS1>;

/**
 * @brief Type for MultipatchType: A field of spline coefficients for a non-batched spline defined
 * on the first of the Patch's logical dimensions.
 */
template <class Patch>
using ConstSplineCoeff1OnPatch_1D = DConstField<typename Patch::IdxRangeBS1>;



// DERIVATIVES -----------------------------------------------------------------------------------

/**
 * @brief Type for MultipatchType: A constant field of the n-th derivatives in the direction of Patch's first
 * logical grid defined on the index range of n and the Patch's second logical grid.
 */
template <class Patch>
using ConstDeriv1_OnPatch_2D
        = DConstField<IdxRange<ddc::Deriv<typename Patch::Dim1>, typename Patch::Grid2>>;

/**
 * @brief Type for MultipatchType: A constant field of the n-th derivatives in the direction of Patch's second
 * logical grid defined on the index range of n and the Patch's first logical grid.
 */
template <class Patch>
using ConstDeriv2_OnPatch_2D
        = DConstField<IdxRange<typename Patch::Grid1, ddc::Deriv<typename Patch::Dim2>>>;

/**
 * @brief Type for MultipatchType: A constant field of the derivatives \partial_1^(n)\partial_2^(m) f(x)
 * defined on the index ranges of valid n and m.
 */
template <class Patch>
using ConstDeriv12_OnPatch_2D
        = DConstField<IdxRange<ddc::Deriv<typename Patch::Dim1>, ddc::Deriv<typename Patch::Dim2>>>;
