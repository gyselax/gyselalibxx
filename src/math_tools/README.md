# Utility Functions

This folder contains mathematical classes and functions.

## Derivative tools

Functions for calculating derivatives with different methods:

- `spline_partial_derivatives.hpp` : calculate derivatives using spline interpolation
- `central_fdm_partial_derivatives.hpp` : calculate derivatives using Finite Difference Method. At the boundary, order one left and right FDM derivatives are used

## Utility tools

The l\_norm\_tools.hpp file contains functions computing the infinity norm. For now, it computes the infinity norm of 
- a double: $`\Vert x \Vert_{\infty} = x`$; 
- a coordinate: $`\Vert x \Vert_{\infty} = \max_{i} (|x_i|)`$.
