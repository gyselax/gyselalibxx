@page interpolation Interpolation Methods

Interpolation methods are any methods which, given values @f$f_j=f(x_j)@f$, allow the construction of an interpolating function @f$\phi(x)@f$ such that:

@f$\phi(x_j)=f_j=f(x_j)@f$

and return the value of that interpolating function @f$\phi(x)@f$ at the desired coordinates.

The general interface for this method is defined via IInterpolator.

The interpolation methods implemented are:
-  SplineInterpolator

## Spline Interpolation

Interpolation by a spline interpolating function is implemented in the class SplineInterpolator. In this case the interpolating function @f$\phi(x)@f$ is a spline. The basis splines are passed as a template parameter, it is these that determine the order of the spline interpolation.

In order for the interpolation to function correctly the values @f$f_j=f(x_j)@f$ provided must be located at the points @f$x_j@f$ identified as the spline interpolation points.

The spline interpolation method is based entirely on the SplineBuilder and SplineEvaluator classes which are found in @ref sll.

## Memory concerns

SplineInterpolator contains a 1D array of spline coefficients. These are unused in most of the code but are used repeatedly in advections. As a result the class PreallocatableSplineInterpolator exists (which inherits from the more general IPreallocatableInterpolator). This class allows a SplineInterpolator to be allocated locally. It is stored in an InterpolatorProxy which means it is deallocated once it goes out of scope. This ensures that the 1D array is not occupying memory during the execution of the rest of the code, but also that the array is only allocated once in each advection operator.
