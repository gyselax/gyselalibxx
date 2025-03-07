# Utility Functions

This folder contains mathematical classes and functions.

## Derivative tools

Functions for calculating derivatives with different methods:

- `spline_1d_partial_derivative.hpp` : calculate derivatives using a 1d spline representation.
- `spline_2d_partial_derivative.hpp` : calculate derivatives using a 2d spline representation.

## Utility tools

The l\_norm\_tools.hpp file contains functions computing the infinity norm. For now, it computes the infinity norm of 
- a double: $`\Vert x \Vert_{\infty} = x`$; 
- a coordinate: $`\Vert x \Vert_{\infty} = \max_{i} (|x_i|)`$.

The spline_builder_2d_cache.hpp file contains a class to be used with a 2d spline representation to compute partial derivatives. The cache class allows for calling the 2d splines builder only once per iteration, even if partial derivatives are evaluated in two directions.
