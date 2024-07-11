# Commented example: the uniform heat equation {#first_steps}
<!--
Copyright (C) The DDC development team, see COPYRIGHT.md file

SPDX-License-Identifier: MIT
-->

In \subpage uniform_heat_equation "examples/uniform_heat_equation.cpp" is a DDC example implementing a forward
finite-difference solver for the heat equation over a rectangle 2D domain with periodic boundary
conditions.

As usual, the file starts with a few includes that will be used in the code.

\snippet uniform_heat_equation.cpp includes

As you can see, to use DDC, we have to include `<ddc/ddc.hpp>`

# Definition of the discretization

Before solving the equation, DDC's primary goal is to define a discrete domain of dimensions specified by the user, along with a discretization along each dimension, which also needs to be specified. 

![domains_image](./images/domains.png "Domains")

Each point in the *DiscreteDomain* is a *DiscreteElement*. These concepts will be clarified later. Let's start by constructing this DiscreteDomain necessary for solving our 2D problem for the heat equation.


## Dimensions naming

We start by defining types that we later use to name our
dimensions.

First, we create `X`: a type that is declared but never defined to act as a tag identifying our
first dimension.

\snippet uniform_heat_equation.cpp X-dimension

Then, we create `DDimX`, a type that will also act as a tag, but we define it to be a uniform
discretization of `X`.
Thus `DDimX` is an identifier for a discrete dimension.

\snippet uniform_heat_equation.cpp X-discretization

We do the same thing for the second dimension `Y`.

\snippet uniform_heat_equation.cpp Y-space

And once again, now for the time dimension.

\snippet uniform_heat_equation.cpp time-space


## Domains

### Dimension X

Once the types are defined, we can start the `main` function where we will define our various
domains. Here for each dimension, the user needs to specify the starting and ending coordinates of the dimension of the domain, as well as the number of discretization points along each of these dimensions. Additionally, we specify here the physical characteristics specific to our equation (the thermal diffusion coefficient).

\snippet uniform_heat_equation.cpp main-start-x-parameters

We start by defining `gwx`, the number of "ghost" points in the `X` dimension, those represented in
dark grey in the figure on each side of the domain.

\snippet uniform_heat_equation.cpp X-parameters
The type for this constant is `DiscreteVector<DDimX>` that represents a number of elements in the
discretization of the `X` dimension.

Once done, we initialize the discretization of the `X` dimension.
Its name has been specified before, but we now set its parameters with the `init_discretization`
function.

\snippet uniform_heat_equation.cpp X-global-domain
Depending on the way the function is called, its return type can differ.
Here we use it with an inner call to `init_ghosted` and receive four 1D domains as a result.
Their type is not specified because we use C++
[structured bindings](https://en.cppreference.com/w/cpp/language/structured_binding), but they are
all of the same type: `DiscreteDomain<DDimX>` that represents a set of contiguous points in the
discretization of `X`.

\ref ddc::UniformPointSampling::init_ghosted "init_ghosted" takes as parameters the coordinate of the first and last discretized points, the
number of discretized points in the domain and the number of additional points on each side of the
domain.
The fours `DiscreteDomain`s returned are:
* `x_domain`: the main domain from the start (included) to end (included) but excluding "ghost"
  points.
* `ghosted_x_domain`: the domain including all "ghost" points.
* `x_pre_ghost`: the "ghost" points that come before the main domain.
* `x_post_ghost`: the "ghost" points that come after the main domain.

The parameters of raw C++ types like `double` or `size_t` can not be used as-is since DDC enforces
strong typing.
Instead, for a coordinate in the `X` dimension, we use `Coordinate<X>` and --as already mentioned--
for a number of elements in the discretization of `X`, we use `DiscreteVector<DDimX>`.

Once this is done, we define two additional domains that will be mirrored to the ghost at the start and at the end of the domain.

\snippet uniform_heat_equation.cpp X-domains

`x_domain_begin` is the sub-domain at the beginning of `x_domain` of the same shape as
`x_post_ghost` that will be mirrored to it.
Reciprocally for `x_domain_end` with `x_pre_ghost`.

The type of both sub-domains is `DiscreteDomain<DDimX>` even if only `DiscreteDomain` is specified,
this relies on C++
[class template argument deduction (CTAD)](https://en.cppreference.com/w/cpp/language/class_template_argument_deduction).
The first parameter given to the constructor is the first element of the domain, a
`DiscreteElement<DDimX>`, the second parameter is a number of elements to include, a
`DiscreteVector<DDimX>`.

**To summarize,** in this section, we have introduced the following types:
* `Coordinate<X>` that represents a point in the continuous dimension `X`,
* `DiscreteElement<DDimX>` that represents one of the elements of `DDimX`, the discretization of
  `X`,
* `DiscreteVector<DDimX>` that represents a number of elements of `DDimX`,
* `DiscreteDomain<DDimX>` that represents an interval in `DDimX`, a set of contiguous elements of
  the  discretization.


### Dimension Y

The domains in the `Y` dimension are handled in a way very similar to the `X` dimension. We first define the domain characteristics along this dimension 

\snippet uniform_heat_equation.cpp main-start-y-parameters

Then we initialize the domain along this dimension just like we did with the `X` dimension.

\snippet uniform_heat_equation.cpp Y-domains

### Time dimension

Then we handle the domains for the simulated time dimension. We first give the simulated time at which to stard and end the simulation. 

\snippet uniform_heat_equation.cpp main-start-t-parameters

Then we use the CFL condition to determine the time step of the simulation.

\snippet uniform_heat_equation.cpp CFL-condition

Finally, we determine the number of time steps and as we did with the `X` and `Y` dimensions, we create the time domain. 

\snippet uniform_heat_equation.cpp time-domain


# Data allocation

We allocate two 2D Chunks along the X and Y dimensions which will be used to map temperature to the domains' points at t and t+dt. When constructing the Chunks one can give an optional string to label the memory allocations. This helps debugging and profiling applications using the Kokkos tools, see also [Kokkos Tools](https://github.com/kokkos/kokkos-tools).

These chunks map the temperature into the full domain (including ghosts) twice:
+ *ghosted_last_temp* for the last fully computed time-step.
+ *ghosted_next_temp* for time-step being computed.

\snippet uniform_heat_equation.cpp data allocation

Note that the `DeviceAllocator` is responsible for allocating memory on the default memory space.

# Initial conditions

To set the initial conditions, the `ghosted_intial_temp` is created and acts as a pointer to the chunk. The const qualifier makes it clear that ghosted_initial_temp always references the same chunk, `ghosted_last_temp` in this case.

\snippet uniform_heat_equation.cpp initial-chunkspan

Then, we iterate over each *DiscreteElement* of the domain to fill `ghosted_initial_temp` with the initial values of the simulation.

\snippet uniform_heat_equation.cpp fill-initial-chunkspan

To display the data, a chunk is created on the host.

\snippet uniform_heat_equation.cpp host-chunk

We deepcopy the data from the `ghosted_last_temp` chunk to `ghosted_temp` on the host.

\snippet uniform_heat_equation.cpp initial-deepcopy

And we display the initial data.

\snippet uniform_heat_equation.cpp initial-display


\snippet uniform_heat_equation.cpp time iteration

To display the data, a chunk is created on the host.


\snippet uniform_heat_equation.cpp host-chunk

We deepcopy the data from the `ghosted_last_temp` chunk to `ghosted_temp` on the host.

\snippet uniform_heat_equation.cpp boundary conditions

\snippet uniform_heat_equation.cpp initial-deepcopy

And we display the initial data.

\snippet uniform_heat_equation.cpp initial-display

For the numerical scheme, two chunkspans are created: 
+ `next_temp` a span excluding ghosts of the temperature at the time-step we will build.
+ `last_temp` a read-only view of the temperature at the previous time-step.Note that *span_cview* returns a read-only ChunkSpan.

\snippet uniform_heat_equation.cpp manipulated views

We then solve the equation.

\snippet uniform_heat_equation.cpp numerical scheme
# Time loop

\snippet uniform_heat_equation.cpp time iteration


## Periodic conditions

\snippet uniform_heat_equation.cpp output

\snippet uniform_heat_equation.cpp boundary conditions


## Numerical scheme


\snippet uniform_heat_equation.cpp swap

For the numerical scheme, two chunkspans are created: 
+ `next_temp` a span excluding ghosts of the temperature at the time-step we will build.
+ `last_temp` a read-only view of the temperature at the previous time-step.Note that *span_cview* returns a read-only ChunkSpan.

\snippet uniform_heat_equation.cpp manipulated views

We then solve the equation.


\snippet uniform_heat_equation.cpp final output

\snippet uniform_heat_equation.cpp numerical scheme
