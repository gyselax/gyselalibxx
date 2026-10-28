# Interpolation Methods

Interpolation methods are any methods which, given values $`f_j=f(x_j)`$, allow the construction of an interpolating function $\phi(x)$ such that:

$`\phi(x_j)=f_j=f(x_j)`$

which can return the value of that interpolating function $\phi(x)$ at the desired coordinates.

This folder describes an operator IInterpolator which takes the values of the function to be interpolated and returns its values at the provided points.

The interpolation methods implemented are:

- SplineInterpolator
- LagrangeInterpolator

The folder also contains tools useful for constructing and evaluating interpolating functions. The sub-folder [polar\_splines](./polar_splines/README.md) contains methods specific to polar splines.

An operator IInterpolator2D is also available for interpolating on a 2D domain.

The interpolation methods implemented are:

- SplineInterpolator2D

## Spline Interpolation

Interpolation by a spline interpolating function is implemented in the class SplineInterpolator. In this case the interpolating function $\phi(x)$ is a spline. The basis splines are passed as a template parameter, it is these that determine the order of the spline interpolation.

In order for the interpolation to function correctly the values $`f_j=f(x_j)`$ provided must be located at the points $`x_j`$ identified as the spline interpolation points.

The spline interpolation method is based entirely on the SplineBuilder and SplineEvaluator classes which are found in DDC.

### Composite Spline Builders

Multi-dimensional splines are tensor products of 1D splines, so their coefficients can be computed by building the spline one dimension at a time. The classes `CompositeSplineBuilder2D` and `CompositeSplineBuilder3D` use this property to combine a 1D builder along the first dimension with a DDC builder for the remaining dimensions:

- `CompositeSplineBuilder2D<Builder1, Builder2>` combines two 1D builders (`Builder2` is a `ddc::SplineBuilder`).
- `CompositeSplineBuilder3D<Builder1, Builder23>` combines a 1D builder with a 2D builder (`Builder23` is a `ddc::SplineBuilder2D`).

The builder along the first dimension only needs to satisfy the `concepts::InterpolationBuilder1D` concept. This makes it possible to use an alternative method in this dimension, e.g. an operator which distributes the computation of the coefficients with MPI.

The composite builders take the same arguments as the equivalent DDC builders (`ddc::SplineBuilder2D` and `ddc::SplineBuilder3D`) and give the same coefficients.

### Polar Spline Interpolation

There is no method to construct a polar spline from the values of a function. It should be possible to construct such a `PolarSplineBuilder`, but it is not clear where the interpolation points should be placed near the O-point in order to obtain a well-conditioned problem. The B-splines, splines and the spline evaluator for the polar splines can be found in the sub-folder [polar\_splines](./polar_splines/README.md).
