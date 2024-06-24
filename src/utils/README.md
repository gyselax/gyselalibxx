# Utility Functions

This folder contains classes and functions which facilitate the writing of the rest of the code. For the most part these objects are required to fill gaps in DDC which will hopefully be filled in the future.

The class ddcHelper exists to provide functionalities which are currently missing from DDC.

The class NDTag exists to provide a way to group directional tags together. This is notably useful in order to create a vector field.

## Utility tools

The utils\_tools.hpp file contains functions computing the infinity norm. For now, it computes the infinity norm of 
- a double: $`\Vert x \Vert_{\infty} = x`$; 
- a coordinate: $`\Vert x \Vert_{\infty} = \max_{i} (|x_i|)`$.
