# Utility Functions

This folder contains mathematical classes and functions.

## Derivative tools

Functions for calculating derivatives with different methods:

- `spline_1d_partial_derivative.hpp` : calculate derivatives using a 1d spline representation.

## Utility tools

The l\_norm\_tools.hpp file contains functions computing the infinity norm. For now, it computes the infinity norm of 
- a double: $`\Vert x \Vert_{\infty} = x`$; 
- a coordinate: $`\Vert x \Vert_{\infty} = \max_{i} (|x_i|)`$.
