# Spline interpolator in polar coordinates

Use a function basis built from a product of two 1D spline bases in $r$ and in $\theta$. This interpolator does **not** use the polar spline basis with special treatment for the center. This means that care must be taken when using this interpolator, to ensure that numerically identical values are provided at $r=0$. If this rule is respected then the representation will be continuous.

The interpolator is built from a product of b-splines if $r$ and $\theta$ rather than using a polar spline basis as the former is significantly less costly than the latter. Therefore the polar spline basis should only be used when additional smoothness is required.


The `operator()` takes as input the values of the function we want to interpolate and the coordinates where we want to interpolate. It computes the coefficients on the 2D bsplines basis and evaluates the decomposition on the coordinates given as input. 
