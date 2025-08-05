

# File types.hpp

[**File List**](files.md) **>** [**data\_types**](dir_2cbcac1ff0802c0a6551cceb4db325f2.md) **>** [**types.hpp**](types_8hpp.md)

[Go to the documentation of this file](types_8hpp.md)


```C++
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

template <class Patch>
using Grid1OnPatch = typename Patch::Grid1;

template <class Patch>
using Grid2OnPatch = typename Patch::Grid2;


// FIELDS ----------------------------------------------------------------------------------------

// --- ON DEVICE ---
template <class Patch>
using DFieldMemOnPatch = DFieldMem<typename Patch::IdxRange12>;

template <class Patch>
using DFieldMem1OnPatch = DFieldMem<typename Patch::IdxRange1>;

template <class Patch>
using DFieldMem2OnPatch = DFieldMem<typename Patch::IdxRange2>;


template <class Patch>
using DFieldOnPatch = DField<typename Patch::IdxRange12>;

template <class Patch>
using DField1OnPatch = DField<typename Patch::IdxRange1>;


template <class Patch>
using DConstFieldOnPatch = DConstField<typename Patch::IdxRange12>;

template <class Patch>
using DConstField1OnPatch = DConstField<typename Patch::IdxRange1>;


// --- ON HOST ---
template <class Patch>
using DFieldOnPatch_host = host_t<DFieldOnPatch<Patch>>;

template <class Patch>
using DConstFieldOnPatch_host = host_t<DConstFieldOnPatch<Patch>>;


// VECTOR FIELDS ---------------------------------------------------------------------------------

// --- ON DEVICE ---
template <class Patch>
using DVectorFieldMemOnPatch = VectorFieldMem<
        double,
        typename Patch::IdxRange12,
        VectorIndexSet<typename Patch::Dim1, typename Patch::Dim2>>;

template <class Patch>
using DVectorFieldOnPatch = VectorField<
        double,
        typename Patch::IdxRange12,
        VectorIndexSet<typename Patch::Dim1, typename Patch::Dim2>>;

template <class Patch>
using DVectorConstFieldOnPatch = VectorConstField<
        double,
        typename Patch::IdxRange12,
        VectorIndexSet<typename Patch::Dim1, typename Patch::Dim2>>;


// IDX, IDXRANGE ---------------------------------------------------------------------------------

template <class Patch>
using IdxRangeOnPatch = typename Patch::IdxRange12;

template <class Patch>
using IdxRange1OnPatch = typename Patch::IdxRange1;

template <class Patch>
using Idx1OnPatch = typename Patch::Idx1;


// IDXRANGESLICE ---------------------------------------------------------------------------------

template <class Patch>
using IdxRange1SliceOnPatch = IdxRangeSlice<typename Patch::Grid1>;

template <class Patch>
using IdxRange2SliceOnPatch = IdxRangeSlice<typename Patch::Grid2>;


// COORDINATE FIELDS -----------------------------------------------------------------------------

// --- ON DEVICE ---
template <class Patch>
using CoordFieldMemOnPatch = FieldMem<typename Patch::Coord12, typename Patch::IdxRange12>;

template <class Patch>
using CoordFieldOnPatch = Field<typename Patch::Coord12, typename Patch::IdxRange12>;

template <class Patch>
using CoordConstFieldOnPatch = ConstField<typename Patch::Coord12, typename Patch::IdxRange12>;


template <class Patch>
using Coord1Field1OnPatch_1D = Field<typename Patch::Coord1, typename Patch::IdxRange1>;


// SPLINES ---------------------------------------------------------------------------------------

template <class Patch>
using BSplines1OnPatch = typename Patch::BSplines1;

template <class Patch>
using BSplines2OnPatch = typename Patch::BSplines2;


// --- ON DEVICE ---
template <class Patch>
using SplineCoeffMemOnPatch_2D = DFieldMem<typename Patch::IdxRangeBS12>;

template <class Patch>
using SplineCoeffOnPatch_2D = DField<typename Patch::IdxRangeBS12>;

template <class Patch>
using ConstSplineCoeffOnPatch_2D = DConstField<typename Patch::IdxRangeBS12>;


template <class Patch>
using SplineCoeff1OnPatch_2D = DField<IdxRange<typename Patch::BSplines1, typename Patch::Grid2>>;

template <class Patch>
using ConstSplineCoeff1OnPatch_2D
        = DConstField<IdxRange<typename Patch::BSplines1, typename Patch::Grid2>>;


template <class Patch>
using SplineCoeff1OnPatch_1D = DField<typename Patch::IdxRangeBS1>;

template <class Patch>
using ConstSplineCoeff1OnPatch_1D = DConstField<typename Patch::IdxRangeBS1>;


// --- ON HOST ---
template <class Patch>
using SplineCoeffMemOnPatch_2D_host = host_t<SplineCoeffMemOnPatch_2D<Patch>>;

template <class Patch>
using SplineCoeffOnPatch_2D_host = host_t<SplineCoeffOnPatch_2D<Patch>>;

template <class Patch>
using ConstSplineCoeffOnPatch_2D_host = host_t<ConstSplineCoeffOnPatch_2D<Patch>>;


// DERIVATIVES -----------------------------------------------------------------------------------

// --- ON DEVICE ---
template <class Patch>
using ConstDeriv1_OnPatch_2D
        = DConstField<IdxRange<ddc::Deriv<typename Patch::Dim1>, typename Patch::Grid2>>;

template <class Patch>
using ConstDeriv2_OnPatch_2D
        = DConstField<IdxRange<typename Patch::Grid1, ddc::Deriv<typename Patch::Dim2>>>;

template <class Patch>
using ConstDeriv12_OnPatch_2D
        = DConstField<IdxRange<ddc::Deriv<typename Patch::Dim1>, ddc::Deriv<typename Patch::Dim2>>>;


// FIELDS WITH DERIVATIVES -----------------------------------------------------------------------

// --- ON DEVICE ---
template <class Patch>
using DerivFieldMemOnPatch = DerivFieldMem<
        double,
        IdxRange<
                ddc::Deriv<typename Patch::Dim1>,
                typename Patch::Grid1,
                ddc::Deriv<typename Patch::Dim2>,
                typename Patch::Grid2>,
        1>;

template <class Patch>
using DerivFieldOnPatch = DerivField<
        double,
        IdxRange<
                ddc::Deriv<typename Patch::Dim1>,
                typename Patch::Grid1,
                ddc::Deriv<typename Patch::Dim2>,
                typename Patch::Grid2>>;


// --- ON HOST ---
template <class Patch>
using DerivFieldMemOnPatch_host = host_t<DerivFieldMemOnPatch<Patch>>;

template <class Patch>
using DerivFieldOnPatch_host = host_t<DerivFieldOnPatch<Patch>>;
```


