# Interpolation on multipatch geometry

This directory defines structures and methods to apply interpolation operators on every patch of a multipatch geometry.
The multipatch interpolation operators implemented are

- [Multipatch interpolation builder](#multipatch-interpolation-builder) - Definition of `InterpolationBuilderCollection` and `InterpolationBuilder2DCollection`.
- [Multipatch interpolation evaluator](#multipatch-interpolation-evaluator) - Definition of `InterpolationEvaluator2DCollection`.

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

`InterpolationBuilderCollection` is used in the same way with 1D builders.

## Multipatch interpolation evaluator

The `InterpolationEvaluator2DCollection` allows the interpolations defined on every patch to be evaluated in one single line.

Any evaluator satisfying `concepts::InterpolationEvaluator` with 2 interpolation dimensions can be used (e.g. a 2D `ddc::SplineEvaluatorND` or `NDLagrangeEvaluator`).
See [interpolation](./../../interpolation/README.md) for more details about the evaluators.
A different evaluator type can be used on each patch.
The class is templated on a `ddc::TypeSeq` of the patches, a `ddc::TypeSeq` of the evaluators (in the same order as the patches), a patch locator and an extrapolation rule.
The evaluators are copied inside the class so that they can be used on GPU.

The patch locator (see [Connectivity](./../connectivity/README.md)) is used to determine on which patch a given coordinate is physically located.
The evaluator of this patch is then used to evaluate the interpolation.
The extrapolation rule determines the value which is returned if the coordinate is outside the domain.
Multipatch extrapolation rules are defined in [spline](./../spline/README.md#multipatch-extrapolation-rules).

The interpolation coefficients are stored in `MultipatchField` objects.
The coordinates are stored in fields stored in `MultipatchField` objects.
Each field corresponds to the coordinates stored on a given patch.
The coordinates are not required to be physically located on the storing patch.

> For example, these coordinates can be characteristic feet stored on the patch where the initial coordinates were.
> Here the coordinates are in the same field if their initial coordinates were on the same patch. They are on the
> same storing patch, but some of them can be physically located outside the storing patch if the feet have crossed
> the edges.

The methods implemented are

- `operator()` to compute a single value or fields of values (stored in `DField` or `DerivField` objects).
- `deriv()` to compute derivatives on the first and/or second dimension on a single coordinate or fields of coordinates.
  No extrapolation rule is applied for derivatives, so the program aborts if a coordinate is outside the domain.
- `integrate()` to compute the integral on each patch. The integrals are stored in a `Kokkos::View` defined on host and of the same size as the number of patches.
  This method is only available if all the evaluators are spline evaluators (`ddc::SplineEvaluatorND`).

**Warning:** The current version of `InterpolationEvaluator2DCollection` does not work on batched domains.

**Warning:** `InterpolationEvaluator2DCollection` only works for patches defined on the same continuous dimensions.

### Example of use of InterpolationEvaluator2DCollection

We define all the evaluators for each patch,

```cpp
EvaluatorType1 evaluator_1(extrapolation_rules_1...);
EvaluatorType2 evaluator_2(extrapolation_rules_2...);
```

We can instantiate an `InterpolationEvaluator2DCollection` from these evaluators, a patch locator and a multipatch extrapolation rule,

```cpp
InterpolationEvaluator2DCollection<
        ddc::TypeSeq<Patch1, Patch2>,
        ddc::TypeSeq<EvaluatorType1, EvaluatorType2>,
        PatchLocator,
        ExtrapolationRule> evaluator(evaluator_1, evaluator_2, patch_locator, extrapolation_rule);
```

We evaluate the interpolations at the coordinates stored on every patch using,

```cpp
evaluator(values, coords, coeffs);
```

## Contents

- `interpolation_builder_2d_collection.hpp`: Define the `InterpolationBuilder2DCollection` operator to apply 2D interpolation builders on every patch.
- `interpolation_builder_collection.hpp`: Define the `InterpolationBuilderCollection` operator to apply 1D interpolation builders on every patch.
- `interpolation_evaluator_2d_collection.hpp`: Define the `InterpolationEvaluator2DCollection` operator to apply 2D interpolation evaluators on every patch.
- `multipatch_interpolation_2d.hpp`: Define the `MultipatchInterpolation2D` class which owns the 2D interpolations of every patch and the matching builder and evaluator collections.
- `multipatch_lagrange_interpolation.hpp`: Define the `MultipatchLagrangeInterpolator2D` alias for a `MultipatchInterpolation2D` using Lagrange interpolation on every patch.
- `multipatch_spline_interpolation.hpp`: Define the `MultipatchSplineInterpolator2D` alias for a `MultipatchInterpolation2D` using spline interpolation on every patch.
