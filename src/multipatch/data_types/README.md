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

Objects which are not plain fields (e.g. `DerivField`, `VectorField`) are stored with
`MultipatchFieldOf` and `MultipatchFieldMemOf`. These are templated on a type template which
gives the type stored on each patch, e.g. `MultipatchFieldOf<DerivFieldOnPatch, Patch1, Patch2>`.

### Types

In `types.hpp` file different aliases are defined to use with `MultipatchType` and `MultipatchFieldOf`:

- To store 1D grids: `Grid1OnPatch` and `Grid2OnPatch`.
- To store fields: `DFieldOnPatch` (and `DFieldOnPatch_host`).
- To store vector fields: `DVectorFieldMemOnPatch`, `DVectorFieldOnPatch` and `DVectorConstFieldOnPatch`.
- To store indices and index ranges:
  - 2D: `IdxRangeOnPatch` and `IdxRangeBSOnPatch` (B-splines);
  - 1D on first dimension: `IdxRange1OnPatch` and `Idx1OnPatch`;
  - Index range slices: `IdxRange1SliceOnPatch` and `IdxRange2SliceOnPatch`.
- To store the spline grids: `BSplines1OnPatch` and `BSplines2OnPatch`.
- To store derivatives:
  - $`\partial_^{(i)} f(x, y_j)`$: `ConstDeriv1_OnPatch_2D`
  - $`\partial_^{(i)} f(x_j, y)`$: `ConstDeriv2_OnPatch_2D`
  - $`\partial_^{(i)} \partial_y^{(j)} f(x, y)`$: `ConstDeriv12_OnPatch_2D`
- To store fields with derivatives: `DerivFieldMemOnPatch`, `DerivFieldOnPatch` and
  `DerivConstFieldOnPatch` (and their `_host` versions).
