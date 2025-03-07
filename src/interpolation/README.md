# Interpolation Methods

Interpolation methods are any methods which, given values $`f_j=f(x_j)`$, allow the construction of an interpolating function $\phi(x)$ such that:

$`\phi(x_j)=f_j=f(x_j)`$

which can return the value of that interpolating function $\phi(x)$ at the desired coordinates.

This folder describes an operator IInterpolator which takes the values of the function to be interpolated and returns its values at the provided points.

The interpolation methods implemented are:
-  SplineInterpolator
-  LagrangeInterpolator

The folder also contains tools useful for constructing and evaluating interpolating functions. The sub-folder [polar\_splines](./polar_splines/README.md) contains methods specific to polar splines.

An operator IInterpolator2D is also available for interpolating on a 2D domain.

The interpolation methods implemented are:
-  SplineInterpolator2D

## Spline Interpolation

Interpolation by a spline interpolating function is implemented in the class SplineInterpolator. In this case the interpolating function $\phi(x)$ is a spline. The basis splines are passed as a template parameter, it is these that determine the order of the spline interpolation.

In order for the interpolation to function correctly the values $`f_j=f(x_j)`$ provided must be located at the points $`x_j`$ identified as the spline interpolation points.

The spline interpolation method is based entirely on the SplineBuilder and SplineEvaluator classes which are found in DDC.

### Polar Spline Interpolation
There is no method to construct a polar spline from the values of a function. It should be possible to construct such a `PolarSplineBuilder`, but it is not clear where the interpolation points should be placed near the O-point in order to obtain a well-conditioned problem. The B-splines, splines and the spline evaluator for the polar splines can be found in the sub-folder [polar\_splines](./polar_splines/README.md).

## Memory concerns

SplineInterpolator contains a 1D array of spline coefficients. These are unused in most of the code but are used repeatedly in advections. As a result the class PreallocatableSplineInterpolator exists (which inherits from the more general IPreallocatableInterpolator). This class allows a SplineInterpolator to be allocated locally. It is stored in an InterpolatorProxy which means it is deallocated once it goes out of scope. This ensures that the 1D array is not occupying memory during the execution of the rest of the code, but also that the array is only allocated once in each advection operator.
