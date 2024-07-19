# Using DDC in Gyselalibxx

[DDC](https://github.com/Maison-de-la-Simulation/ddc) is a library which aims to provide types which represent mathematical/physical concepts.
Representing these concepts with types allows the compiler to enforce the mathematical validity of expressions.

The DDC library is based on templates. The template parameters are based on physical dimensions.

When using DDC the first step is therefore to create structures representing each of the physical dimensions, e.g. RDimX. These objects should only contain one attribute, a static constexpr boolean called PERIODIC which indicates whether the dimension is periodic or not.
The DDC types are then parametrised by this structure. The "R" in the name of the dimension stands for "real", it indicates that the dimension is in real space and is continuous.

The following sections describe some of the DDC types used in Gyselalibxx.

## Coordinate

A `ddc::Coordinate` is one of the only DDC types which represents a continuous data type. This means that it can take any value that can be represented by a double. It represents the position of a coordinate in the vector space.

Coordinates can have 1 or more dimension. E.g. the coordinate of a position on a $(r,\theta)$ slice should have the type `ddc::Coordinate<RDimR, RDimP>`, where `RDimR` represents the radial dimension $r$, and `RDimP` represents the poloidal dimension $\theta$.

If the value of the coordinate needs to be used in a mathematical expression, the scalar (`double`) quantity stored in one of the dimensions of a coordinate can be extracted using `ddc::get<RDimOfInterest>(my_coord)`.

It is also possible to extract a coordinate on a subset of the original dimensions using `ddc::select<RDimOfInterest>(my_coord)`. For example if we want to get the position of an object on a radial slice $(r,\theta)$, but we are given the coordinate in the full vector space $`(r, \theta, \varphi, v_\parallel, \mu)`$ then we can do:
```cpp
ddc::Coordinate<RDimR, RDimP, RDimT, RDimV, RDimM> full_coord(...);
ddc::Coordinate<RDimR, RDimP> slice_coord = ddc::select<RDimR, RDimP>(full_coord);
```

Coordinates can be combined using operators. For example, let us consider three vectors $P$, $Q$, and $R$ defined on a cartesian space $(x,y)$:

![Vector Image](./images/Coordinate_operations.jpg)

The vector $Q$ can be written as $Q=R-P$. Similarly in the code, we would have:
```cpp
ddc::Coordinate<RDimX, RDimY> P(5.0, 2.0);
ddc::Coordinate<RDimX, RDimY> R(7.0, 6.0);
ddc::Coordinate<RDimX, RDimY> Q = R - P;
```

The `Coordinate` class also provides an addition operator, comparison operators, and an output operator for easy printing using `cout`.

In Gyselalibxx the alias `CoordX` is usually defined in `geometry.hpp` to describe the type `ddc::Coordinate<RDimX>` more succinctly.

## Domain

The physical problems that our simulations describe are defined on a domain. The domain on which the problem is defined is continuous (e.g. a radial domain $[0,1)$). However a simulation evolves on a discrete domain. This means that the value of the function is only known at a discrete set of points. In the case of a function $f(x)\rightarrow y \in \mathbb{R}$ with $x \in [0,1)$, we would usually discretise the domain $[0,1)$ as follows:
$$
x_0, x_1, ..., x_i, ..., x_N
$$

DDC provides multiple types to represent the concepts required to interact with such a domain.

### PointSampling

The points $`\{x_0, ..., x_N\}`$ are a point sampling. This sampling can either be uniform or non-uniform. Accordingly DDC provides the 2 classes:
- UniformPointSampling
- NonUniformPointSampling

A uniform point sampling is a collection of points which are equidistant, it is therefore defined with an origin and a step. In contrast the points found in a non-uniform point sampling are arbitrary. This kind of sampling must therefore be initialised from a list of points.

In Gyselalibxx the alias `IDimX` is usually defined in `geometry.hpp` to describe a point sampling along the dimension X. The "I" in `IDimX` stands for "Interpolation". This is because the function is only stored at the interpolation points. The point sampling contains these interpolation points.

### DiscreteElement

The simplest type to understand is `ddc::DiscreteElement`. This corresponds to the index of the point in the point sampling. E.g. the point `x_i` in the point sampling can be indexed using the object `ddc::DiscreteElement<IDimX>(i)`.
A `ddc::DiscreteElement` is therefore roughly equivalent to an integer. Compared to an integer, it additionally contains information about the physical direction being examined. This allows the compiler to raise errors if typos/copy-paste errors lead to the wrong dimension being used.

As a discrete concept this class does not take the continuous dimension as a template parameter, rather it takes the discrete point sampling instead.

In Gyselalibxx the alias `IndexX` is usually defined in `geometry.hpp` to describe the index of a point sampling along the dimension X.

We can also create multi-dimensional indices. E.g. the point `(x_i, y_j)` can be indexed using the object `ddc::DiscreteElement<IDimX, IDimY>(i, j)`

It is also possible to extract an index on a subset of the original dimensions using `ddc::select<IDimOfInterest>(my_nd_index)`.

If a point sampling has been initialised using the function `ddc::init_discrete_space`, a `ddc::DiscreteElement` defined on that sampling can be used to find the coordinates of the points in the point sampling. This is done using the function `ddc::coordinate`.

E.g:
```cpp
ddc::init_discrete_space<ddc::UniformPointSampling<RDimR>>(0.0, 0.1);
ddc::DiscreteElement<IDimR> i(0); // The first point in the point sampling
std::cout << ddc::coordinate(i) << std::endl;
```
will output:
```
(0.0)
```

### DiscreteVector

A `ddc::DiscreteVector` describes the number of grid points between points in a sampling. E.g. `x_3` and `x_6` are separated by `ddc::DiscreteVector<IDimX>(3)`.

As a discrete concept this class does not take the continuous dimension as a template parameter, rather it takes the discrete point sampling instead.

This type is useful when we have the index of a point and we need to get the next point. A `ddc::DiscreteVector` can be added or subtracted from a `ddc::DiscreteElement` as long as both objects are templated by the same dimension. Similarly a `ddc::DiscreteVector` is the result of subtracting 2 `ddc::DiscreteElement`s.

E.g:
```cpp
ddc::DiscreteElement<IDimX> i(4);
ddc::DiscreteElement<IDimX> j(6);
ddc::DiscreteVector<IDimX> k = j-i;
i += k;
```

In Gyselalibxx the alias `IVectX` is usually defined in `geometry.hpp` to describe the vector from one element of a point sampling along the dimension X to another.

As with `ddc::DiscreteElement`s, a `ddc::DiscreteVector` can be multi-dimensional and lower dimension `ddc::DiscreteVector` objects can be extracted using `ddc::select<IDimOfInterest>(my_nd_vector)`.

### DiscreteDomain

The last concept necessary to define a discrete domain is the concept of sub-domains. The class `ddc::DiscreteDomain` is designed to describe a sub-domain (although it can of course hold the whole domain). It is templated over each of the discrete point samplings used to describe the discrete domain on which the simulation evolves.

Each subdomain is described by:
-   An origin : This is the `ddc::DiscreteElement` which indicates the first point in the domain.
-   A size : This is a `ddc::DiscreteVector` indicating the number of elements in each dimension.

For example if we consider the 2D domain: $`[x_0, ..., x_N] \times [y_0, ... y_M]`$, the domain would be described as:
```cpp
ddc::DiscreteElement<IDimX, IDimY> origin(0, 0);
ddc::DiscreteVector<IDimX, IDimY> size(N, M);
ddc::DiscreteDomain<IDimX, IDimY> domain(origin, size);
```

Similarly the sub-domain: $`[x_i, ..., x_j] \times [y_k, ... y_l]`$ would be described as:
```cpp
ddc::DiscreteElement<IDimX> i_index(i);
ddc::DiscreteElement<IDimX> j_index(j);
ddc::DiscreteElement<IDimY> k_index(k);
ddc::DiscreteElement<IDimY> l_index(l);
ddc::DiscreteElement<IDimX, IDimY> origin(i_index, k_index);
ddc::DiscreteVector<IDimX, IDimY> size((j_index-i_index), (l_index-k_index));
ddc::DiscreteDomain<IDimX, IDimY> domain(origin, size);
```

When working with domains we do not usually know if we have access to all of a domain or simply a subdomain. It is therefore important to use the domain functions to traverse the domain rather than initialising elements manually as we don't know the index of the first element of a `ddc::DiscreteDomain` at compile time.

There are multiple functions available for traversing a domain. Most of the time we will traverse the entire domain. This can be done simply as `ddc::DiscreteDomain` implements the functions `begin()` and `end()`. These functions are called automatically using the modern C++ syntax for a for element in list, or using the `ddc::for_each` function. The latter is to be preferred as it will allow us to add parallelism later. The syntax is:
```cpp
for (ddc::DiscreteElement<IDimX> index : domain) {
}
```
or:
```cpp
ddc::for_each(domain, [&](ddc::DiscreteElement<IDimX> index) {
});
```
In the case of a `ddc::for_each` the second argument is a lambda function. The `[&]` ensures that any variable defined outside the loop are captured by reference so they can be used inside the lamda function.

It is also common to need to iterate over a sub-domain. Sub-domains can be created using the syntax described above, however `ddc::DiscreteDomain` also contains several helper functions which are designed to facilitate the creation of subdomains:
- `take_first(ddc::DiscreteVector<..> n)` : Returns a sub-domain, restricted to the first n elements.
- `take_last(ddc::DiscreteVector<..> n)` : Returns a sub-domain, restricted to the last n elements.
- `remove_first(ddc::DiscreteVector<..> n)` : Returns a sub-domain, containing all elements except the first n elements.
- `remove_last(ddc::DiscreteVector<..> n)` : Returns a sub-domain, containing all elements except the last n elements.
- `remove(ddc::DiscreteVector<..> n_first, ddc::DiscreteVector<..> n_last)` : Returns a sub-domain, containing all elements except the first n\_first elements and the last n\_last elements.

It is also possible to extract a sub-domain on a subset of the original dimensions using `ddc::select<IDimOfInterest>(my_nd_domain)`.

Finally it may not be possible to express the elements you want to iterate over as a domain. This is notably the case if you want to iterate over every j-th element. In this case it is necessary to fall back on `ddc::DiscreteDomainIterator`. The syntax in this case is:
```cpp
for (ddc::DiscreteDomainIterator<IDimX> it=domain.begin(); it < domain.end(); it += j) {
    ddc::DiscreteElement<IDimX> index = *it;
}
```

In addition to the iteration functionalities, `ddc::DiscreteDomain` also has other useful functions. The following is a non-exhaustive list of useful functions:
- `front()` : Returns the first `ddc::DiscreteElement` in the domain.
- `back()` : Returns the last `ddc::DiscreteElement` in the domain.
- `size()` : Returns the total number of points in the domain (the product of the number of points in each dimension).
- `extents()` : Returns the number of points in each dimension stored in a  `ddc::DiscreteVector`.

In Gyselalibxx the alias `IDomainX` is usually defined in `geometry.hpp` to describe the domain or sub-domain containing points from the point sampling `IDimX`.

## Data Storage

DDC provides the type `ddc::Chunk` to store data. This type is parametrised by the underlying data type (e.g. `double`), and the `ddc::DiscreteDomain` on which the values are defined.

In order to initialise the data storage to the correct size, a `ddc::Chunk` is initailised by providing the `ddc::DiscreteDomain` on which the values are defined.

In Gyselalibxx the alias `FieldX` is usually defined in `geometry.hpp` to describe a function defined on the domain or sub-domain containing points from the point sampling `IDimX`. E.g:
```cpp
FieldXVx<double> distribution_function_2d(dom_x_vx);
```

In Gyselalibxx the functions are almost all defined using real numbers. The additional alias `DFieldX` is therefore defined to represent a real function defined on the domain or sub-domain containing points from the point sampling `IDimX`. E.g:
```cpp
DFieldRPTVM distribution_function_2d(dom_radial_poloidal_toroidal_velocity_mu);
```

A `ddc::Chunk` is indexed using either multiple 1D `ddc::DiscreteElement`s, or one ND `ddc::DiscreteElement`. As the `ddc::DiscreteElement`s contain information about the dimension on which they act it is not necessary to pass these arguments in a specific order. Thus the following two functions are equivalent:
```cpp
double get_element_1(DFieldXVx my_chunk, IndexX i, IndexY j) {
    return my_chunk(i, j);
}

double get_element_2(DFieldXVx my_chunk, IndexX i, IndexY j) {
    return my_chunk(j, i);
}
```
This is particularly useful if we don't know the layout order of the data and will allow us to reorder this data without changing the way we interact with the chunk.

Copying a `ddc::Chunk` is potentially costly. In order to avoid accidental copying DDC is structured such that the only way to copy data from one `ddc::Chunk` to another is using the function `ddc::deepcopy`. Using an assignment or initialising from another `ddc::Chunk` will result in a compiler error.

To avoid copying data unnecessarily, DDC provides the type `ddc::ChunkSpan`. This can be thought of as a reference to a `ddc::Chunk`. However it is slightly more complex than this as a `ddc::ChunkSpan` does not have to reference the entire domain stored in the `ddc::Chunk`. The `[]` operator can be passed a `ddc::DiscreteDomain` to create a `ddc::ChunkSpan` which only references part of the `ddc::Chunk`. This is especially useful for accessing slices. In this case a simpler syntax exists, where we only need to pass the index of the slice. For example, if we wish to get a $(r, \theta)$ slice from a distribution function defined in $`(r, \theta, \varphi, v_\parallel, \mu)`$ we would do the following:
```cpp
ddc::ChunkSpan<double, ddc::DiscreteDomain<IDimR, IDimP>> get_slice(ddc::Chunk<double, ddc::DiscreteDomain<IDimR, IDimP, IDimT, IDimV, IDimM>>& distribution_function, ddc::DiscreteElement<IDimT, IDimV, IDimM> index)
{
    return distribution_function[index];
}
```

Unless you need to allocate data, you should always use `ddc::ChunkSpan` rather than `ddc::Chunk`.

In Gyselalibxx the alias `SpanX` is usually defined in `geometry.hpp` to describe a reference to a function defined on the domain or sub-domain containing points from the point sampling `IDimX`. E.g:
```cpp
FieldXVx<double> distribution_function_2d(dom_x_vx);
SpanXVx<double> distribution_function_2d_ref(distribution_function_2d);
```

As for `ddc::Chunk` an additional alias `DSpanX` is defined to simplify the notation for doubles.

Finally we also define the aliases `ViewX` and `DViewX` to represent constant versions of `SpanX` and `DSpanX`. These objects are used for function arguments when the contents of the array must not be altered. As they are aliases, not new objects, all the same functions that work for `ddc::ChunkSpan` can also be used for `ddc::ChunkView`. For example, if we wish to get a $(r, \theta)$ constant slice from a distribution function defined in $`(r, \theta, \varphi, v_\parallel, \mu)`$ we would do the following:
```cpp
ddc::ChunkView<double, ddc::DiscreteDomain<IDimR, IDimP>> get_slice(ddc::Chunk<double, ddc::DiscreteDomain<IDimR, IDimP, IDimT, IDimV, IDimM>>& distribution_function, ddc::DiscreteElement<IDimT, IDimV, IDimM> index)
{
    return distribution_function[index];
}
```

## Example

Let us consider the following system of equations (a simple Vlasov-Poisson system):
$$
\left[\partial_t + \frac{1}{\sqrt{m_s}} \left(v\partial_x - q_s \partial_x \phi(t,x)\partial_v\right)\right] f_s(t,x,v) = 0\\
\partial_x^2\phi(t,x) = - \sum_s q_s \int f_s(t,x,v) dv\\
f_s(0,x,v) = \frac{n_0}{2\pi T_0}\exp\left(-\frac{v^2}{T_0}\right)\cos(k_xx)
$$

which evolves over a 2 dimensional domain defined by the dimensions $(x,v)$. The domain is $[0, 700]\times[-6,6]$.
We can also consider the species s as a discrete dimension which is defined as [i, e].

The first thing which is necessary to define our system are types describing the continuous dimensions:
```cpp
// Tag to represent the x-dimension
struct RDimX {
    static bool constexpr PERIODIC = true;
}

// Tag to represent the v-dimension
struct RDimVx {
    static bool constexpr PERIODIC = false;
}
```

We also need types to define the discrete domain on which the simulation will evolve. The distribution function $`f_s(t,x,v)`$ will evolve over the discrete domain $`[i, e]\times[x_0,...,x_N]\times[v_0,...,v_{N_v}]`$. Point samplings are required to define the positions of the grid points in each of the three dimensions:
- The object $`[x_0,...,x_N]`$ is defined with a point sampling and will be denoted $IDimX$.
- The object $`[v_0,...,v_{N_v}]`$ is defined with a point sampling and will be denoted $IDimVx$.
- The object $['i', 'e']$ is defined as a `SpeciesInformation` collection and will be denoted $IDimSp$.

With these objects defined, the domain(s) can then be created. The domain for the distribution function has the type: `ddc::DiscreteDomain<IDimSp, IDimX, IDimVx>` however for simplicity it is denoted `DomainSpXVx`. The electric potential $\phi$ has a smaller domain with the type `ddc::DiscreteDomain<IDimX>` denoted `DomainX`.

Other domains may be useful to work on slices of data or on subdomains (e.g. if we distribute our function using MPI).

The data itself will be stored in `ddc::Chunk`s. For example the distribution function will be stored in a `ddc::Chunk<double, ddc::DiscreteDomain<IDimSp, IDimX, IDimVx>>`. To improve readability this type is denoted `DFieldSpXVx`. Similarly the electric potential will be stored in a `DFieldX`.

The other types will be useful when interacting with the data. For example, in order to initialise the distribution function we need to loop over the domain using the coordinates and a `ddc::DiscreteElement<IDimSp, IDimX, IDimVx>` (denoted `IndexSpXVx` for simplicity) to fill the array:
```cpp
DomainSpXVx domain(..);
DFieldSpXVx distribution_function(domain);

ddc::for_each(domain, [&](IndexSpXVx index) {
    CoordXVx coord = ddc::coordinate(ddc::select<IDimX, IDimVx>(index));
    double v_pos = ddc::get<RDimVx>(coord);
    double x_pos = ddc::get<RDimX>(coord);
    distribution_function(index) = n0 / (2 * M_PI * T0) * exp(-v_pos * v_pos / T0) * cos(kx * x_pos);
});
```

`ddc::ChunkSpan` objects are required for function arguments. For example the prototype of the function which defines the electric potential would be:
```cpp
DSpanX operator()(DSpanX electric_potential, DViewSpXVx distribution_function);
```
where `DSpanX` is a reference to the `DFieldX` object where the electric potential is stored, and `DViewSpXVx` is a constant reference to the `DFieldSpXVx` object where the distribution function is stored.

### Continuous vs. discrete objects

The majority of the DDC types represent objects in or on a discrete vector space. It may not be immediately obvious why continous objects are also needed given that simulations are necessarily discrete. Backward semi-Lagrangian advection is a good example of where this is required. In backward semi-Lagrangian advection we work on a discrete grid (a `ddc::DiscreteDomain`), for each point of the grid we trace the characteristic to find where a particle positioned at this grid point would have been located at the previous time step. This location, known as the foot of the characteristic can be located anywhere in the vector space, unlike the grid points, the foot of the characteristic is not restricted to a subset of valid positions. It must therefore be represented with a continuous object. In this case a `ddc::Coordinate` perfectly describes the position of the foot of the characteristic while also providing all pertinent information about the dimensions.

### Example use of DiscreteVector for Finite Differences

The class `ddc::DiscreteVector` is particularly useful in cases where we wish to express a stencil. Let us consider the case of a 2D finite differences scheme on uniform meshes expressed as follows:
$$
\nabla^2 f(x,y) = \frac{f(x-h_x, y) + f(x+h_x, y) + f(x, y-h_y) + f(x, y+h_y) - 4f(x,y)}{h_x h_y}
$$

This code can be written simply using `ddc::DiscreteElement` and `ddc::DiscreteVector`:
```cpp
double get_laplacian_at_position(DFieldXY function_values, IndexXY position)
{
    VectX x_step(1);
    VectY y_step(1);

    // Get the uniform point sampling in the appropriate direction and use it to
    // extract the step h.
    double h_x = ddc::discrete_space<RDimX>().step();
    double h_y = ddc::discrete_space<RDimY>().step();

    return (function_values(position - x_step) + function_values(position + x_step)
            + function_values(position - y_step) function_values(position + y_step)
            - 4 * function_values(position))
           / (h_x * h_y);
}
```

As you can see DDC takes care of combining elements and vectors to ensure that the resulting `IndexXY` is in the correct position. This makes it simpler to write such stencils and easier to spot errors, especially in cases involving lots of dimensions.
