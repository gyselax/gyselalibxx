# Data Types for Multipatch Geometry

This directory contains classes which are useful for handling objects and
types defined on a multipatch domain.

## MultipatchType

The class `MultipatchType` stores different objects of a type which is templated with
patches (see `Patch`). The `MultipatchType::get` method can be used to
retrieve an object defined on a specified patch.
For example we can use index ranges on different patches which would be the type
`Patch::IdxRange12`.

So after defining

```cpp
template<class Patch>
using IdxFieldOnPatch = typename Patch::IdxRange12;
```

we could then have three index ranges `idx_range1`, `idx_range2` and `idx_range3` on
patches 1,2 and 3 respectively. The `MultipatchType` object would then
be initialised as

```cpp
MultipatchType<IdxFieldOnPatch, Patch1, Patch2, Patch3> multipatch_idx_range(idx_range1, idx_range2, idx_range3);
```

and the index range on patch 3 can be retrieved via

```cpp
typename Patch3::IdxRange12 index_range3_from_multipatch = multipatch_idx_range.get<Patch3>();
```

## MultipatchField and MultipatchFieldMem

Fields defined on each patch are stored in a `MultipatchField` (references to the data) or a
`MultipatchFieldMem` (which owns the memory). These aliases are templated in the same way as
`Field` and `FieldMem`, with the index range replaced by a `MultipatchType` of index ranges:

```cpp
using MultipatchIdxRange = MultipatchType<IdxRangeOnPatch, Patch1, Patch2>;

MultipatchIdxRange idx_ranges(idx_range_1, idx_range_2);
DMultipatchFieldMem<MultipatchIdxRange> fields_alloc(idx_ranges);
DMultipatchField<MultipatchIdxRange> fields(fields_alloc);
DMultipatchConstField<MultipatchIdxRange> const_fields(fields.get_const_field());
MultipatchConstField<Coord<R, Theta>, MultipatchIdxRange> coords(coords_1, coords_2);
```

Vector fields whose components are the same on every patch are stored in a `MultipatchVectorField`
or a `MultipatchVectorFieldMem`. These are templated in the same way as `VectorField` and
`VectorFieldMem`:

```cpp
DMultipatchVectorFieldMem<MultipatchIdxRange, VectorIndexSet<R, Theta>> vec_fields_alloc(idx_ranges);
DMultipatchVectorField<MultipatchIdxRange, VectorIndexSet<R, Theta>> vec_fields(vec_fields_alloc);
DMultipatchField<MultipatchIdxRange> vec_fields_r = ddcHelper::get<R>(vec_fields);
```

Fields with derivatives are stored in a `MultipatchDerivField` or a `MultipatchDerivFieldMem`.
These are templated in the same way as `DerivField` and `DerivFieldMem`. As for `DerivField`, the
index ranges include the derivative dimensions (see `IdxRangeDerivOnPatch`). The memory can be
allocated from the physical index ranges, in which case the derivatives are defined on the
boundaries of each patch:

```cpp
using MultipatchIdxRangeDeriv = MultipatchType<IdxRangeDerivOnPatch, Patch1, Patch2>;

DMultipatchDerivFieldMem<MultipatchIdxRangeDeriv, 1> functions_and_derivs_alloc(idx_ranges);
DMultipatchDerivField<MultipatchIdxRangeDeriv> functions_and_derivs(functions_and_derivs_alloc);
```

Other objects (e.g. vector fields whose components differ from one patch to another) are stored
with `MultipatchFieldOf` and `MultipatchFieldMemOf`. These are templated on a type template which
gives the type stored on each patch, e.g. `MultipatchFieldOf<DerivFieldOnPatch, Patch1, Patch2>`.
They are also needed where a type template is expected, e.g. the `ValuesOnPatch` template
parameter of `MultipatchSplineEvaluator2D`.

All these aliases describe the same `detail::MultipatchField` and `detail::MultipatchFieldMem`
classes. A function which needs to know the patches can therefore deduce them from a
`detail::MultipatchField<ddc::TypeSeq<Patches...>, TypeSeqInternalTypes>` parameter.

### Types

In `types.hpp` file different aliases are defined to use with `MultipatchType` and `MultipatchFieldOf`:

- To store 1D grids: `Grid1OnPatch` and `Grid2OnPatch`.
- To store fields: `DFieldOnPatch` (and `DFieldOnPatch_host`).
- To store indices and index ranges:
  - 2D: `IdxRangeOnPatch` and `IdxRangeBSOnPatch` (B-splines);
  - 1D on first dimension: `IdxRange1OnPatch` and `Idx1OnPatch`;
  - Index range slices: `IdxRange1SliceOnPatch` and `IdxRange2SliceOnPatch`;
  - Index range of a field with derivatives, on (dx, x, dy, y): `IdxRangeDerivOnPatch`;
  - Index ranges of derivatives: `IdxRangeDeriv1Grid2OnPatch` on (dx, y), `IdxRangeGrid1Deriv2OnPatch`
    on (x, dy) and `IdxRangeDeriv12OnPatch` on (dx, dy).
- To store the spline grids: `BSplines1OnPatch` and `BSplines2OnPatch`.
- To store fields with derivatives: `DerivFieldOnPatch`.
