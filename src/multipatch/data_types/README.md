# Data Types for Multipatch Geometry

This directory contains classes which are useful for handling objects and
types defined on a multipatch domain.

## MultipatchType

The class `MultipatchType` stores different objects of a type which is templated with
patches (see `Patch`). The `MultipatchType::get` method can be used to
retrieve an object defined on a specified patch.
For example we can use fields on different patches which would be the type
`DField<Patch::IdxRange12>`.

So after defining

```
template<class Patch>
using DFieldOnPatch = DField<Patch::IdxRange12>;
```

we could then have three fields `field1`, `field2` and `field3` on
patches 1,2 and 3 respectively. The `MultipatchType` object would then
be initialised as

```
MultipatchType<DFieldOnPatch, Patch1, Patch2, Patch3> multipatch_field(field1, field2, field3);
```

and the field on patch 3 can be retrieved via

```
DField<Patch3::IdxRange12> field3_from_multipatch = multipatch_field.get<Patch3>();
```

### Types

In `types.hpp` file different aliases are defined to use the MultipatchType class:

* To store 1D grids: `Grid1OnPatch` and `Grid2OnPatch`.
* To store fields and constant fields:
  * 2D: (`DFieldMemOnPatch`), `DFieldOnPatch` and `DConstFieldOnPatch`;
  * 1D on first dimension: `DField1OnPatch` and `DConstField1OnPatch`.
* To store indices and index ranges:
  * 2D: `IdxRangeOnPatch`
  * 1D on first dimension: `IdxRange1OnPatch` and `Idx1OnPatch`.
* To store coordinates:
  * 2D: `CoordFieldOnPatch`
  * 1D on first dimension: `Coord1Field1OnPatch_1D`.
* To store the splines:
  * Spline grids: `BSplines1OnPatch` and `BSplines2OnPatch`.
  * 2D: 2D spline coefficients `SplineCoeffOnPatch_2D` and `ConstSplineCoeffOnPatch_2D`;
  * 2D: spline coefficients on the first dimension `SplineCoeff1OnPatch_2D` and `ConstSplineCoeff1OnPatch_2D`;
  * 1D on the fist dimension: spline coefficients `SplineCoeff1OnPatch_1D` and `SplineCoeff2OnPatch_1D`;
* To store derivatives:
  * $`\partial_^{(i)} f(x, y_j)`$: `ConstDeriv1_OnPatch_2D`
  * $`\partial_^{(i)} f(x_j, y)`$: `ConstDeriv2_OnPatch_2D`
  * $`\partial_^{(i)} \partial_y^{(j)} f(x, y)`$: `ConstDeriv12_OnPatch_2D`
