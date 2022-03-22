# Commented example: the heat equation {#first_steps}

In \subpage heat_equation "examples/heat_equation.cpp" is a DDC example implementing a forward
finite-difference solver for the heat equation over a rectangle 2D domain with periodic boundary
conditions.

As usual, the file starts with a few includes that will be used in the code.

\snippet heat_equation.cpp includes

As you can see, DDC includes all follow the same convention: `<ddc/SYMBOL>` where `SYMBOL` is a
the name of a DDC symbol.
So for example, in order to use a class named `Chunk`, you should include `<ddc/Chunk>` and to use a
free function template named `for_each`, you should include `<ddc/for_each>`.

Then, we define the value of some parameters that would typically be read from some form of
configuration file in a more realistic code.

\snippet heat_equation.cpp parameters


# Definition of the discretization

Given the previous parameters, we will now define a 2D point discretization.

![domains_image](./images/domains.png "Domains")


## Dimensions naming

Before starting the `main` function, we start by defining types that we later use to name our
dimensions.

First, we create `X`: a type that is declared but never defined to act as a tag identifying our
first dimension.

\snippet heat_equation.cpp X-dimension

Then, we create `DDimX`, a type that will also act as a tag, but we define it to be a uniform
discretization of `X`.
Thus `DDimX` is an identifier for a discrete dimension.

\snippet heat_equation.cpp X-discretization

We do the same thing for the second dimension `Y`.

\snippet heat_equation.cpp Y-space

And once again, now for the time dimension.

\snippet heat_equation.cpp time-space


## Domains

### Dimension X

Once the types are defined, we can start the `main` function where we will define our various
domains.

\snippet heat_equation.cpp main-start

We start by defining `gwx`, the number of "ghost" points in the `X` dimension.
This can be made `static constexpr` since here this is a compile-time constant.
The type for this constant is `DiscreteVector<DDimX>` that represents a number of elements in the
discretization of the `X` dimension.

\snippet heat_equation.cpp X-parameters

Once done, we initialize the actual discretization of the `X` dimension.
Its name has been specified before, but we now set its parameters.

\snippet heat_equation.cpp X-global-domain

Since we use the `init_ghosted` function, we receive four 1D domains as a result.
Their type is not specified because we use C++
[structured bindings](https://en.cppreference.com/w/cpp/language/structured_binding), but they are
all of the same type: `DiscreteDomain<DDimX>` that represents a set of contiguous points in the
discretization of `X`.
This function expects as parameter the coordinate of the first and last discretized points, the
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

Once this is done, we define two additional domains.

\snippet heat_equation.cpp X-domains

`x_domain_begin` is the sub-domain at the beginning of `x_domain` of the same shape as
`x_post_ghost` that will be mirrored to it.
Reciprocally for `x_domain_end` with `x_pre_ghost`.

The type of both sub-domains is `DiscreteDomain<DDimX>` even if only `DiscreteDomain` is specified,
this relies on C++
[class template argument deduction (CTAD)](https://en.cppreference.com/w/cpp/language/class_template_argument_deduction).
The first parameter given to the constructor is the first element of the domain, a
`DiscreteCoordinate<DDimX>`, the second parameter is a number of elements to include, a
`DiscreteVector<DDimX>`.

**To summarize,** in this section, we have introduced the following types:
* `Coordinate<X>` that represents a point in the continuous dimension `X`,
* `DiscreteCoordinate<DDimX>` that represents one of the elements of `DDimX`, the discretization of
  `X`,
* `DiscreteVector<DDimX>` that represents a number of elements of `DDimX`,
* `DiscreteDomain<DDimX>` that represents an interval in `DDimX`, a set of contiguous elements of
  the  discretization.


### Dimension Y

The domains in the `Y` dimension are handled in a way very similar to the `X` dimension.

\snippet heat_equation.cpp Y-domains


### Time dimension

Then we handle the domains for the simulated time dimension.

\snippet heat_equation.cpp time-domains


# Data allocation

We allocate two 2D Chunks along the X and Y dimensions which will be used to map temperature to the domains' points at t and t+dt.

\snippet heat_equation.cpp data allocation


# Initial conditions

\snippet heat_equation.cpp initial-conditions

\snippet heat_equation.cpp initial output


# Time loop

\snippet heat_equation.cpp time iteration


## Periodic conditions

\snippet heat_equation.cpp boundary conditions


## Numerical scheme

\snippet heat_equation.cpp manipulated views

\snippet heat_equation.cpp numerical scheme


## Output

\snippet heat_equation.cpp output


## Final swap

\snippet heat_equation.cpp swap


# Final output

\snippet heat_equation.cpp final output
