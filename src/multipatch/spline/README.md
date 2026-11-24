# Spline on multipatch geometry

This directory defines structures and methods to deal with different splines on every patches.
The multipatch patch operators implemented are

- [Multipatch extrapolation rules](#multipatch-extrapolation-rules) - Definition of extrapolation rules depending on the geometries for `InterpolationEvaluator2DCollection`.

The spline builders and evaluators can be applied on every patch using `InterpolationBuilderCollection`, `InterpolationBuilder2DCollection` and `InterpolationEvaluator2DCollection` (see [multipatch interpolation](./../interpolation/README.md)).

## Multipatch extrapolation rules

To evaluate at a given coordinate outside the domain, extrapolation rules have been implemented. They are declined according to the geometry:

- `NullExtrapolationRule`: general for every geometry and every interpolation method. It sets the values at zero for every coordinates outside the domain.

- `ConstantExtrapolationRuleOnion`: specialised for onion geometries (see `OnionPatchLocator` and [Connectivity](./../connectivity/README.md)) and spline interpolations. It calls `ddc::ConstantExtrapolationRule` to evaluate the outside coordinate.
 Two areas are considered as outside: with radius inferior to the minimum radius and with radius superior to the maximum radius.

## Contents

- `constant_extrapolation_rules_onion.hpp`: Define `ConstantExtrapolationRuleOnion`, constant extrapolation rule for onion geometries.
- `null_extrapolation_rules.hpp`: Define `NullExtrapolationRule`, a null extrapolation rule for every geometry.  
