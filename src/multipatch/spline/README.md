# Spline on multipatch geometry

This directory defines structures and methods to deal with different splines on every patches.
The multipatch patch operators implemented are 

* [Multipatch spline builder](#src_multipatch_spline__Multipatch_spline_builder) - Definition of `MultipatchSplineBuilder` and `MultipatchSplineBuilder2D`. 
* [Multipatch spline evaluator](#src_multipatch_spline__Multipatch_spline_evaluator) - Definition of `MultipatchSplineEvaluator2D`. 
* [Multipatch extrapolation rules](#src_multipatch_spline__Multipatch_extrapolation_rules) - Definition of extrapolation rules depending on the geometries for `MultipatchSplineEvaluator2D`. 


## Multipatch spline builder 

The `MultipatchSplineBuilder` and `MultipatchSplineBuilder2D` allow all the given spline builders to be called in one single line. 

The template parameters of this class reproduce the `SplineBuilder` and `SplineBuilder2D` classes from DDC. 
The spline coefficients and values of the functions are stored in `MultipatchType` objects. 
The spline builders are stored inside the class in a `std::tuple`.

For Hermite boundary conditions, derivatives have to be provided. They are also stored in `MultipatchType` objects. 

### Example of use of MultipatchSplineBuilder
We define all the builders for each patch, 
```cpp
SplineBuilderType1 builder_1(idx_range_1); 
SplineBuilderType2 builder_2(idx_range_2); 
SplineBuilderType3 builder_3(idx_range_3); 
...
```

We can instantiate a MultipatchSplineBuilder from these builders, 
```cpp
MultipatchSplineBuilder builder (builder_1, builder_2, builder_3); 
```

We allocate memory for the spline representations on each patch,
```cpp
SplineType1 spline_1(spline_idx_range_1); 
SplineType2 spline_2(spline_idx_range_2); 
SplineType3 spline_3(spline_idx_range_3); 
...
MultipatchType<SplineTypeOnPatch, Patch1, Patch2, Patch3> patches_splines (spline_1, spline_2, spline_3); 
```

We do the same with the values of the function, 
```cpp
ValuesType1 values_1(idx_range_1); 
ValuesType2 values_2(idx_range_2); 
ValuesType3 values_3(idx_range_3); 
...
// initialisation 
...
MultipatchType<ValuesTypeOnPatch, Patch1, Patch2, Patch3>  patches_values (values_1, values_2, values_3); 
```

And we build the spline representation on every patch using, 
```cpp
builder(patches_splines, patches_values); 
```

## Multipatch spline evaluator

The `MultipatchSplineEvaluator2D` allows all the spline evaluators to be called in one single line. 

Contrary to the `MultipatchSplineBuilder` class, the spline evaluators are not given as input but we directly 
compute the evaluation inside the class. The template parameters and aliases are similar to the ones in the `SplineEvaluator2D`
class in DDC. 

To instantiate a `MultipatchSplineEvaluator2D`, we need a patch locator operator to be able to determine on which patch the given (see [Connectivity](./../connectivity/README.md))
coordinates will be physically located. We also need extrapolation rules to know what type of value returning if case of the 
coordinates are outside of the domain. Multipatch extrapolation rules are implemented and defined below [Multipatch extrapolation rules](#src_multipatch_spline__Multipatch_extrapolation_rules).

The coordinates are stored in fields stored in `MultipatchType` objects.
They are stored on patches. Each field corresponds to the coordinates stored on a given patch. 
The coordinates are not enforced to be physically located on the storing patch.

> For example, these coordinates can be characteristic feet stored on the patch where the initial coordinates were. 
> Here the coordinates are in the same field if their initial coordinates were on the same patch. They are on the 
> same storing patch, but some of them can be physically located outside of the storing patch if the feet have crossed 
> the edges. 


Similarly to `SplineEvaluator2D`, methods to get the derivatives are implemented. No extrapolation rules are given, so `MultipatchSplineEvaluator2D` will throw an exception error if we give a coordinate outside of the domain. 
The methods implemented are
* `operator()` to compute a single value or fields of values. 
* `deriv_dim_1()`, `deriv_dim_2()` to compute derivatives on the first or second dimension on a single coordinates or fields of coordinates. 
* `deriv<InterestDim>()` to compute derivatives on the first or second dimension on a single coordinates only.
* `deriv_1_and_2()` to compute cross-derivatives n a single coordinates or fields of coordinates. 
* `integrate()` to compute the integral on each patch. The integral are stored in a `Kokkos::View` defined on host and of the same size as the number of patches. 


**Warning:** The current version of `MultipatchSplineEvaluator2D` does not work on batched domain. 

**Warning:** The mappings applied in the given patch locator have to contain `operator()` from the logical domain to the physical domain and from the physical domain to the logical domain. Both operators are called in the `MultipatchSplineEvaluator2D` class to compute equivalent coordinates from one patch to another. 


## Multipatch extrapolation rules

To evaluate at a given coordinate outside of the domain, extrapolation rules have been implemented. They are declined according to the geometry: 

 * `NullExtrapolationRule`: general for every geometries. It sets the values at zero for every coordinates outside of the domain. 

 * `ConstantExtrapolationRuleOnion`: specialised for onion geometries (see `OnionPatchLocator` and [Connectivity](./../connectivity/README.md)). It calls `ddc::ConstantExtrapolationRule` to evaluate the outside coordinate. 
 Two areas are considered as outside: with radius inferior to the minimum radius and with radius superior to the maximum radius. 


## Contents
- `constant_extrapolation_rules_onion.hpp`: Define `ConstantExtrapolationRuleOnion`, constant extrapolation rule for onion geometries. 
- `multipatch_spline_builder_2d.hpp`: Define the `MultipatchSplineBuilder2D` operator to apply 2D spline builders on every patches. 
- `multipatch_spline_builder.hpp`: Define the `MultipatchSplineBuilder` operator to apply 1D spline builders on every patches. 
- `multipatch_spline_evaluator_2d.hpp`: Define the `MultipatchSplineEvaluator2D` operator to apply 2D spline evaluators on every patches. 
- `null_extrapolation_rules.hpp`: Define `NullExtrapolationRule`, a null extrapolation rule for every geometries.  




