# Interpolation on multipatch geometry

This directory defines structures and methods to apply interpolation operators on every patch of a multipatch geometry.
The multipatch interpolation operators implemented are

- [Multipatch interpolation builder](#multipatch-interpolation-builder) - Definition of `InterpolationBuilderCollection` and `InterpolationBuilder2DCollection`.

## Multipatch interpolation builder

The `InterpolationBuilderCollection` and `InterpolationBuilder2DCollection` allow all the given 1D or 2D interpolation builders to be called in one single line.

For `InterpolationBuilderCollection`, any builder satisfying `concepts::InterpolationBuilder1D` can be used (e.g. `ddc::SplineBuilder` or `IdentityInterpolationBuilder`).
Each builder must interpolate along one of the grids (`Grid1` or `Grid2`) of its patch. The values may contain batch dimensions.
For `InterpolationBuilder2DCollection`, any builder satisfying `concepts::InterpolationBuilder` with 2 interpolation dimensions can be used (e.g. `ddc::SplineBuilder2D` or `NDIdentityInterpolationBuilder`).
See [interpolation](./../../interpolation/README.md) for more details about the builders.
A different builder type can be used on each patch (e.g. to use different boundary conditions).
The class is templated on a `ddc::TypeSeq` of the patches and a `ddc::TypeSeq` of the builders (in the same order as the patches).
The interpolation coefficients and values of the functions are stored in `MultipatchField` objects.
The builders are stored inside the class in a `std::tuple`.

Derivatives can be provided as optional `MultipatchField` objects. For `InterpolationBuilder2DCollection` they can also be provided by passing the values in `DerivField` objects.
For `InterpolationBuilderCollection`, the derivatives at the lower and upper bounds are passed to every builder (as required by `concepts::InterpolationBuilder1D`).
For `InterpolationBuilder2DCollection`, they are only used by `ddc::SplineBuilder2D` builders (e.g. for Hermite boundary conditions) and are ignored by other builders.

### Example of use of InterpolationBuilder2DCollection

We define all the builders for each patch,

```cpp
BuilderType1 builder_1(idx_range_1);
BuilderType2 builder_2(idx_range_2);
```

We can instantiate an `InterpolationBuilder2DCollection` from these builders,

```cpp
InterpolationBuilder2DCollection<
        ddc::TypeSeq<Patch1, Patch2>,
        ddc::TypeSeq<BuilderType1, BuilderType2>> builder(builder_1, builder_2);
```

We allocate memory for the interpolation coefficients on each patch,

```cpp
MultipatchType<IdxRangeOnPatch, Patch1, Patch2> idx_ranges(idx_range_1, idx_range_2);
MultipatchFieldMem<CoeffsMemOnPatch, Patch1, Patch2> coeffs_alloc(builder.batched_basis_idx_ranges(idx_ranges));
MultipatchField<CoeffsOnPatch, Patch1, Patch2> coeffs(coeffs_alloc);
```

We store the values of the function on each patch,

```cpp
MultipatchField<DConstFieldOnPatch, Patch1, Patch2> values(get_const_field(values_1), get_const_field(values_2));
```

And we build the interpolation coefficients on every patch using,

```cpp
builder(coeffs, values);
```

`InterpolationBuilderCollection` is used in the same way with 1D builders. It does not provide `batched_basis_idx_ranges`, so the coefficients are allocated using `batched_basis_idx_range(builder_i, idx_range_i)` on each patch.

## Contents

- `interpolation_builder_2d_collection.hpp`: Define the `InterpolationBuilder2DCollection` operator to apply 2D interpolation builders on every patch.
- `interpolation_builder_collection.hpp`: Define the `InterpolationBuilderCollection` operator to apply 1D interpolation builders on every patch.
