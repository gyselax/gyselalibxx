// SPDX-License-Identifier: MIT
#pragma once

/*
    Geometry defined here: 9 * 2D patches.
        - for all patches: 
            - non-periodic on X,Y;
            - non-uniform cubic splines;
            - non-uniform Grid1 and Grid2; 
        - across patches:
            - periodic on X;
            - non-periodic on Y;


      1  |  2  |  3
    -----------------
      4  |  5  |  6
    -----------------
      7  |  8  |  9
*/

#define GEOM_NAMESPACE_NAME periodic_strips_non_uniform_2d_9patches
#include "9patches_2d_periodic_strips.hpp"
#undef GEOM_NAMESPACE_NAME
