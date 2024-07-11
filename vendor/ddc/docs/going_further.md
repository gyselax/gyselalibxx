# Commented example: the non uniform heat equation {#going_further}
<!--
Copyright (C) The DDC development team, see COPYRIGHT.md file

SPDX-License-Identifier: MIT
-->

In \subpage non_uniform_heat_equation "examples/non_uniform_heat_equation.cpp" is a DDC example implementing a forward
finite-difference solver for the heat equation over a rectangle 2D domain with periodic boundary
conditions and non-uniform space discretization.

As usual, the file starts with a few includes that will be used in the code.

\snippet non_uniform_heat_equation.cpp includes

# Differences with the uniform problem resolution

For non-uniform discretization, differences arise compared to uniform discretization, particularly in terms of dimension labeling, domain construction, and even in the resolution of the problem.

## Dimensions naming

Just like the uniform case, we start by defining types that we later use to name our
dimensions. 

\snippet non_uniform_heat_equation.cpp X-dimension

However, we then define the discretization as non uniform.

\snippet non_uniform_heat_equation.cpp X-discretization

We do the same thing for the second dimension `Y`.

\snippet non_uniform_heat_equation.cpp Y-space

And for this case, we kept an uniform discretization for the time dimension.

\snippet non_uniform_heat_equation.cpp time-space

## Domains

### Dimension X

Just like the uniform case, we define the main characteristics of the domain. 

\snippet non_uniform_heat_equation.cpp main-start-x-parameters

The first step to create a `DiscreteDomain` with non-uniform spatial discretization is to create a C++ iterator containing the actual coordinates of the points along each dimension. For tutorial purposes, we have constructed a function `generate_random_vector` that generates a vector with random data ranging from the lower bound to the higher one to populate our vectors.

For the `X` dimension, we want to build 4 domains: 
* `x_domain`: the main domain from the start (included) to end (included) but excluding "ghost"
  points.
* `ghosted_x_domain`: the domain including all "ghost" points.
* `x_pre_ghost`: the "ghost" points that come before the main domain.
* `x_post_ghost`: the "ghost" points that come after the main domain.

To do so, we have to create the iterator for the main domain.

\snippet non_uniform_heat_equation.cpp iterator_main-domain

For the initialization of ghost points in the non-uniform case, it was necessary to explicitly describe the position of the pre-ghost and post-ghost points. For the pre-ghost, it was necessary to start from the first coordinate along x and subtract the difference between the last and the penultimate coordinate of the x dimension to maintain periodicity. For the post-ghost point, take the last point and add the difference between the first and second points along the x dimension.

\snippet non_uniform_heat_equation.cpp ghost_points_x

Then we can create our 4 domains using the `init_discretization`
function with the `init_ghosted` function that takes the vectors as parameters. 

\snippet non_uniform_heat_equation.cpp build-domains

**To summarize,** in this section we saw how to build a domain with non-uniform space discretization: 
+ Create an iterator for your domain.
+ Call the `init_discretization` function with `init_ghosted` if you need ghost points or with `init` for a regular domain. 

### Dimension Y 

For the `Y` dimension, we do the same. We start by defining our 3 vectors. 

\snippet non_uniform_heat_equation.cpp Y-vectors 

And we use them to build the 4 domains in the same way as for the X dimension.

\snippet non_uniform_heat_equation.cpp build-Y-domain

### Time dimension

Then we handle the domains for the simulated time dimension. We first give the simulated time at which to stard and end the simulation. 

\snippet non_uniform_heat_equation.cpp main-start-t-parameters

The CFL conditions are more challenging to achieve for the case of non-uniform discretization. 
Since the spatial steps are not uniform, we first need to find the maximum of the inverse of the square of the spatial step for the `X` and `Y` dimensions. And then we obtain the minimum value for `dt`.

\snippet non_uniform_heat_equation.cpp CFL-condition

We can calculate the number of time steps and build the `DiscreteDomain` for the time dimension.

\snippet non_uniform_heat_equation.cpp time-domain

## Time loop

Allocation and initialization are the same as for the uniform case. Let's focus on resolving the numerical scheme.
The main difference in solving the numerical equation is that we need to account for the fact that the values of dx and dy on the left and right sides are different. We use the functions `distance_at_left` and `distance_at_right` to solve the equation. 

\snippet non_uniform_heat_equation.cpp numerical scheme




















