# Quadrature Methods

Quadrature methods are methods for calculating the integral of a function from its value at a given set of points $p_1, ..., p_n$.
These methods can be written as a scalar product of the value of the function at the points $f(p_1), ..., f(p_n)$, with the quadrature coefficients $q_1, ..., q_n$.

This folder provides the class Quadrature which integrates a function from its values at the points defined by the discrete domain(s) on which the class is defined.
The class should be initialised with the quadrature coefficients.

Helper functions provide the quadrature coefficients obtained using different quadrature methods.
The methods currently implemented are:
-  trapezoid_quadrature_coefficients()
-  spline_quadrature_coefficients()


Additionally the function quadrature_coeffs_nd() helps define multi-dimensional quadrature methods from 1D methods.

In the compute_norms.hpp file, the functions compute_L1_norm() and compute_L2_norm() return the norms of a function with a given quadrature method. 
The function compute_coeffs_on_mapping() add the Jacobian determinant of the mapping as factor of the quadrature coefficients. 