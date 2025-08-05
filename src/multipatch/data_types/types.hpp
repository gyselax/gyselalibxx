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
    List of aliases useful to store data in MultipatchType or MultipatchField. 
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
/// @brief A FieldMem of doubles defined on the Patch's 2D logical domain.
template <class Patch>
using DFieldMemOnPatch = DFieldMem<typename Patch::IdxRange12>;

/// @brief A FieldMem of doubles defined on the Patch's first logical domain.
template <class Patch>
using DFieldMem1OnPatch = DFieldMem<typename Patch::IdxRange1>;

/// @brief A FieldMem of doubles defined on the Patch's second logical domain.
template <class Patch>
using DFieldMem2OnPatch = DFieldMem<typename Patch::IdxRange2>;


/// @brief A Field of doubles defined on the Patch's 2D logical domain.
template <class Patch>
using DFieldOnPatch = DField<typename Patch::IdxRange12>;

/// @brief A Field of doubles defined on the Patch's first logical domain.
template <class Patch>
using DField1OnPatch = DField<typename Patch::IdxRange1>;


/// @brief A constant Field of doubles defined on the Patch's 2D logical domain.
template <class Patch>
using DConstFieldOnPatch = DConstField<typename Patch::IdxRange12>;

/// @brief A constant Field of doubles on the Patch's first logical domain.
template <class Patch>
using DConstField1OnPatch = DConstField<typename Patch::IdxRange1>;


// --- ON HOST ---
/// @brief A Field of doubles defined on host on the Patch's 2D logical domain.
template <class Patch>
using DFieldOnPatch_host = host_t<DFieldOnPatch<Patch>>;

/// @brief A constant Field of doubles on host on the Patch's 2D logical domain.
template <class Patch>
using DConstFieldOnPatch_host = host_t<DConstFieldOnPatch<Patch>>;


// VECTOR FIELDS ---------------------------------------------------------------------------------

// --- ON DEVICE ---
/// @brief A VectorFieldMem of doubles defined on the Patch's 2D logical domain.
template <class Patch>
using DVectorFieldMemOnPatch = VectorFieldMem<
        double,
        typename Patch::IdxRange12,
        VectorIndexSet<typename Patch::Dim1, typename Patch::Dim2>>;

/// @brief A VectorField of doubles defined on the Patch's 2D logical domain.
template <class Patch>
using DVectorFieldOnPatch = VectorField<
        double,
        typename Patch::IdxRange12,
        VectorIndexSet<typename Patch::Dim1, typename Patch::Dim2>>;

/// @brief A ConstVectorField of doubles defined on the Patch's 2D logical domain.
template <class Patch>
using DVectorConstFieldOnPatch = VectorConstField<
        double,
        typename Patch::IdxRange12,
        VectorIndexSet<typename Patch::Dim1, typename Patch::Dim2>>;


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


// IDXRANGESLICE ---------------------------------------------------------------------------------

/// @brief An index range slice over the grids on the first Patch's 1D logical domain.
template <class Patch>
using IdxRange1SliceOnPatch = IdxRangeSlice<typename Patch::Grid1>;

/// @brief An index range slice over the grids on the second Patch's 1D logical domain.
template <class Patch>
using IdxRange2SliceOnPatch = IdxRangeSlice<typename Patch::Grid2>;


// COORDINATE FIELDS -----------------------------------------------------------------------------

// --- ON DEVICE ---
/// @brief FieldMem of 2D coordinates defined on the 2D Patch's logical domain.
template <class Patch>
using CoordFieldMemOnPatch = FieldMem<typename Patch::Coord12, typename Patch::IdxRange12>;

/// @brief Field of 2D coordinates defined on the 2D Patch's logical domain.
template <class Patch>
using CoordFieldOnPatch = Field<typename Patch::Coord12, typename Patch::IdxRange12>;

/// @brief ConstField of 2D coordinates defined on the 2D Patch's logical domain.
template <class Patch>
using CoordConstFieldOnPatch = ConstField<typename Patch::Coord12, typename Patch::IdxRange12>;


/// @brief Field of 1D first coordinates defined on the first Patch's 1D logical domain.
template <class Patch>
using Coord1Field1OnPatch_1D = Field<typename Patch::Coord1, typename Patch::IdxRange1>;


// SPLINES ---------------------------------------------------------------------------------------

/// @brief The BSplines dimension on the first of the Patch's logical dimensions.
template <class Patch>
using BSplines1OnPatch = typename Patch::BSplines1;

/// @brief The BSplines dimension on the second of the Patch's logical dimensions.
template <class Patch>
using BSplines2OnPatch = typename Patch::BSplines2;


// --- ON DEVICE ---
/**
 * @brief A field memory of 2D spline coefficients for a non-batched spline defined
 * on both of the Patch's logical dimensions.
 */
template <class Patch>
using SplineCoeffMemOnPatch_2D = DFieldMem<typename Patch::IdxRangeBS12>;

/**
 * @brief A field of 2D spline coefficients for a non-batched spline defined
 * on both of the Patch's logical dimensions.
 */
template <class Patch>
using SplineCoeffOnPatch_2D = DField<typename Patch::IdxRangeBS12>;

/**
 * @brief A constant field of 2D spline coefficients for a non-batched spline defined
 * on both of the Patch's logical dimensions.
 */
template <class Patch>
using ConstSplineCoeffOnPatch_2D = DConstField<typename Patch::IdxRangeBS12>;


/**
 * @brief A field of spline coefficients batched over the second of the Patch's
 * logical dimensions for a spline defined on the first of the Patch's logical dimensions.
 */
template <class Patch>
using SplineCoeff1OnPatch_2D = DField<IdxRange<typename Patch::BSplines1, typename Patch::Grid2>>;

/**
 * @brief A constant field of spline coefficients batched over the second of the Patch's
 * logical dimensions for a spline defined on the first of the Patch's logical dimensions.
 */
template <class Patch>
using ConstSplineCoeff1OnPatch_2D
        = DConstField<IdxRange<typename Patch::BSplines1, typename Patch::Grid2>>;


/**
 * @brief A field of spline coefficients for a non-batched spline defined
 * on the first of the Patch's logical dimensions.
 */
template <class Patch>
using SplineCoeff1OnPatch_1D = DField<typename Patch::IdxRangeBS1>;

/**
 * @brief A constant field of spline coefficients for a non-batched spline defined
 * on the first of the Patch's logical dimensions.
 */
template <class Patch>
using ConstSplineCoeff1OnPatch_1D = DConstField<typename Patch::IdxRangeBS1>;


// --- ON HOST ---
/**
 * @brief A field memory of 2D spline coefficients for a non-batched spline defined
 * on both of the Patch's logical dimensions. Defined on host. 
 */
template <class Patch>
using SplineCoeffMemOnPatch_2D_host = host_t<SplineCoeffMemOnPatch_2D<Patch>>;

/**
 * @brief A field of 2D spline coefficients for a non-batched spline defined
 * on both of the Patch's logical dimensions. Defined on host.
 */
template <class Patch>
using SplineCoeffOnPatch_2D_host = host_t<SplineCoeffOnPatch_2D<Patch>>;

/**
 * @brief A constant field of 2D spline coefficients for a non-batched spline defined
 * on both of the Patch's logical dimensions. Defined on host. 
 */
template <class Patch>
using ConstSplineCoeffOnPatch_2D_host = host_t<ConstSplineCoeffOnPatch_2D<Patch>>;


// FIELDS WITH DERIVATIVES -----------------------------------------------------------------------

// --- ON DEVICE ---
/**
 * @brief A field memory with derivatives constaining the function values on the 2D grid, the first derivatives 
 * on the boundaries of the 2D grid and the cross-derivatives on the corner of the 2D grid for the 2D grid 
 * defined on the given Patch. The data are stored on an index range on (dx, x, dy, y) for a grid on (x,y). 
 */
template <class Patch>
using DerivFieldMemOnPatch = DerivFieldMem<
        double,
        IdxRange<
                ddc::Deriv<typename Patch::Dim1>,
                typename Patch::Grid1,
                ddc::Deriv<typename Patch::Dim2>,
                typename Patch::Grid2>,
        1>;

/**
 * @brief A field with derivatives constaining the function values on the 2D grid, and first derivatives 
 * on the boundaries of the 2D grid and the cross-derivatives on the corner of the 2D grid for the 2D grid 
 * defined on the given Patch. The data are stored on an index range on (dx, x, dy, y) for a grid on (x,y). 
 */
template <class Patch>
using DerivFieldOnPatch = DerivField<
        double,
        IdxRange<
                ddc::Deriv<typename Patch::Dim1>,
                typename Patch::Grid1,
                ddc::Deriv<typename Patch::Dim2>,
                typename Patch::Grid2>>;


// --- ON HOST ---
/**
 * @brief A field memory with derivatives constaining the function values on the 2D grid, and first derivatives 
 * on the boundaries of the 2D grid and the cross-derivatives on the corner of the 2D grid for the 2D grid 
 * defined on the given Patch. The data are stored on an index range on (dx, x, dy, y) for a grid on (x,y). 
 * Defined on host. 
 */
template <class Patch>
using DerivFieldMemOnPatch_host = host_t<DerivFieldMemOnPatch<Patch>>;

/**
 * @brief A field with derivatives constaining the function values on the 2D grid, and first derivatives 
 * on the boundaries of the 2D grid and the cross-derivatives on the corner of the 2D grid for the 2D grid 
 * defined on the given Patch. The data are stored on an index range on (dx, x, dy, y) for a grid on (x,y). 
 * Defined on host. 
 */
template <class Patch>
using DerivFieldOnPatch_host = host_t<DerivFieldOnPatch<Patch>>;
