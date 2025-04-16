# 1D-1V Landau Damping Simulation with Gyselalib++

> **Tutorial Goal:** Learn to use Gyselalib++ to implement a 1D-1V Vlasov–Poisson simulation demonstrating Landau damping, by walking through the scientific and algorithmic choices made in [this simulation](../../simulations/geometryXVx/landau/landau_fft.cpp).

---

## Introduction

In this tutorial we will walk you through the implementation of a **1D-1V Landau damping simulation**.

Each decision we make in building the simulation from the domain geometry to the choice of numerical solver corresponds to a statement you might find in a scientific article describing such a simulation. We show how to express these decisions clearly using Gyselalib++.

---

## Problem Setup: Vlasov–Poisson in 1D-1V

We consider the **electrostatic Vlasov–Poisson system** in one spatial and one velocity dimension. This models the time evolution of the **electron distribution function** $ f(x, v, t) $, which obeys the Vlasov equation:

```math
\frac{\partial f}{\partial t} + v \frac{\partial f}{\partial x} + E(x, t) \frac{\partial f}{\partial v} = 0
```

Here, $ E(x, t) $ is the self-consistent electric field, determined by Poisson’s equation:

```math
\frac{\partial E}{\partial x} = \int_{-\infty}^{\infty} q_s f(x, v, t)\, dv - n_{0}
```

This system conserves phase-space density and exhibits **Landau damping** — the decay of the electric field due to phase mixing in velocity space, even in the absence of collisions.

In Gyselalib++, each of these equations corresponds to a **numerical operator** that we must explicitly construct and compose.

---

## Continuous Geometry

Before discretizing the system, we must define the **computational domain** for both space and velocity.

We consider a **periodic spatial domain**:
- $ x \in [0, L_x] $, with $ L_x = 2\pi $
- Periodic boundary conditions in $ x $

And a **bounded velocity domain**:
- $ v \in [-v_{\max}, v_{\max}] $, with $ v_{\max} = 6 $
- Homogeneous Dirichlet (zero) boundary conditions in $ v $

These boundary choices are typical for Landau damping simulations:
- The periodicity in $ x $ allows us to use **Fourier-based Poisson solvers**
- The truncation of the $ v $-domain must be large enough to contain the significant support of the distribution function

We can now take the first steps to construct our simulation. We create structures to describe the continuous dimensions of the **computational domain**:

```cpp
/**
 * @brief A class which describes the real space in the spatial X direction
 */
struct X
{
    /// @brief A boolean indicating if the dimension is periodic.
    static bool constexpr PERIODIC = true;

    /// A boolean indicating if dimension describes a covariant coordinate.
    static bool constexpr IS_COVARIANT = true;
    /// A boolean indicating if dimension describes a contravariant coordinate.
    static bool constexpr IS_CONTRAVARIANT = true;
    /// A type-alias mapping to the co/contra-variant counterpart.
    using Dual = X;
};

/**
 * @brief A class which describes the real space in the X-velocity direction
 */
struct Vx
{
    /// @brief A boolean indicating if the dimension is periodic.
    static bool constexpr PERIODIC = false;
};
```

---

## Numerical Methods

Our numerical scheme is based on a **splitting approach** to solve the Vlasov–Poisson system. Each term in the equations is discretized and handled by an explicit operator, which mirrors the structure of a simulation methods section.

The key choices are:
- A **semi-Lagrangian method** for phase-space advection
- A **Strang splitting** technique for time integration
- An **FFT-based Poisson solver** for computing the electric field
- A **spline-based quadrature rule** for evaluating the charge density

Each of these is a deliberate, transparent modeling decision — in Gyselalib++, they correspond to clearly defined building blocks in the simulation.

---

### Predictor–Corrector Scheme

To improve the accuracy of the time integration, we use a **predictor–corrector scheme** to solve the system of equations. This involves:
1. Calculating the electric field
2. Computing an intermediate solution using the calculated electric field
3. Calculating the electric field at the intermediate solution
4. Computing the final solution using the electric field calculated at the intermediate solution

This enhances stability and temporal accuracy without the need for small time steps.

---

### Poisson Solver

To compute the electric field, we solve Poisson’s equation:

```math
\frac{d^2 \phi}{dx^2} = \rho(x) - n_0, \quad E(x) = -\frac{d\phi}{dx}
```

We use a **spectral solver based on FFT**, leveraging the periodic boundary condition in space. This provides high efficiency and spectral accuracy for the potential and field.

---

### Charge Density via Quadrature

To compute the charge density $ \rho(x) $, we integrate the distribution function over velocity space:

```math
\rho(x) = \int_{-\infty}^{\infty} f(x, v)\, dv
```

This integral is evaluated using a **high-order quadrature rule** based on a spline approximation.

---

### Operator Splitting

We use **Strang splitting** to decompose the Vlasov equation:

```math
\frac{\partial f}{\partial t} = \mathcal{A}_x(f) + \mathcal{A}_v(f)
```

Where:
- $` \mathcal{A}_x(f) = -v \frac{\partial f}{\partial x} `$
- $` \mathcal{A}_v(f) = -E(x) \frac{\partial f}{\partial v} `$

Strang splitting alternates between these operators in a **symmetric sequence** to achieve second-order accuracy in time:

```math
f^{n+1} = e^{\frac{\Delta t}{2} \mathcal{A}_v} \, e^{\Delta t \mathcal{A}_x} \, e^{\frac{\Delta t}{2} \mathcal{A}_v} f^n
```

This approach is standard in semi-Lagrangian Vlasov solvers, and allows us to treat each advection step independently.

---

### Semi-Lagrangian Advection

Each sub-step of the splitting is handled via a **semi-Lagrangian method**. This involves tracing characteristics **backward in time** to find where information originates, and interpolating from that location.
In our simulation we will use a **spline representation** of $ f $ for the interpolation.

Benefits of this approach:
- Unconditionally stable for large time steps
- Works naturally with non-uniform grids and smooth basis functions
- Avoids CFL restrictions typical of Eulerian methods

When tracing characteristics backward in time the foot of the characteristic may fall outside the simulation domain. We define **extrapolation conditions** to describe how the function is evaluated in these cases:
- In the **spatial domain**, periodic extrapolation is used, consistent with the assumed periodicity of the domain.
- In the **velocity domain**, constant extrapolation is used. As the distribution function decays rapidly at large velocities we can assume that the function is zero outside of the domain, but the use of a constant extrapolation instead avoids the accidental introduction of discontinuities.

### Spline Representation

The semi-Lagrangian method employed in this simulation requires the representation of the distribution function $ f(x, v, t) $ as a one-dimensional function during each substep of the Strang splitting scheme. Specifically, during advection in configuration space, $ f(\cdot, v) $ is treated as a function of $ x $ for fixed $ v $, and vice versa for velocity space.

We represent these one-dimensional functions using **B-spline interpolation**. The spline spaces used in each dimension are defined as follows:

- **Uniform knot spacing** is used in both the spatial and velocity domains.
- **Cubic B-splines** (polynomial degree $ p = 3 $) are employed to provide smooth and accurate interpolation.

This representation is used consistently throughout the simulation for both initialization and interpolation during advection steps.

The spline $S(x)$ is evaluated from the B-splines using the following equation:
```math
S(x) = \sum_i c_i b_i(x)
```
where $`c_i`$ are the spline coefficients and $`b_i(x)`$ are the B-splines.

Building a spline representation of a function is equivalent to calculating the spline coefficients $`c_i`$.
Our simulation will evolve on grid points $`x_i`$ so we can define $`N_p`$ equations:
```math
S(x_j) = \sum_i c_i b_i(x_j)
```
where $`N_p`$ is the number of interpolation points.
For our simulation we choose to use **Greville abscissae** as the interpolation points, as this ensures good approximation properties. However this may result in fewer equations than the degrees of freedom of the spline space. We therefore also define **closure conditions** which reflect the physical properties of each domain:
- In the **spatial domain**, periodic boundary conditions are enforced, consistent with the assumed periodicity of the Vlasov–Poisson system.
- In the **velocity domain**, the spline space enforces **homogeneous Hermite-type boundary conditions**, which ensure that derivatives vanish at the domain boundaries. This reflects the rapid decay of the distribution function at large velocities.

---

## Constructing the Simulation

In this section, we build the Landau damping simulation step by step, beginning with the **spline basis** and **grid**, and working our way up to the full operator composition.

This construction reflects the mathematical structure of the Vlasov–Poisson system and the numerical methods we introduced earlier, from foundational representation to time integration.

Fundamental operators which are useful for many simulations and geometries are provided by Gyselalib++, however as we move towards more specific operators the dependencies on the geometry or the simulation become stronger.
In these cases you may need to write your own operator built upon the tools that Gyselalib++ provides. Nevertheless, as we will see, Gyselalib++ probides several geometry specific operators for common use cases.

---

### Defining the Spline Basis

The spline basis as defined [above](#spline-representation) is described using DDC:

```cpp
int constexpr BSDegreeX = 3;
int constexpr BSDegreeVx = 3;

struct BSplinesX : ddc::UniformBSplines<X, BSDegreeX>
{
};
struct BSplinesVx : ddc::UniformBSplines<Vx, BSDegreeVx>
{
};
```

In the body of the simulation these objects must be initialised as global objects. This is done as follows:
```cpp
Coord<X> const x_min(0.0);
Coord<X> const x_max(1.0);
IdxStep<GridX> const x_ncells(100);

Coord<Vx> const vx_min(-6);
Coord<Vx> const vx_max(6);
IdxStep<GridVx> const vx_ncells(30);

ddc::init_discrete_space<BSplinesX>(x_min, x_max, x_ncells);
ddc::init_discrete_space<BSplinesVx>(vx_min, vx_max, vx_ncells);
```

The function `ddc::init_discrete_space` calls the [constructor of `ddc::UniformBSplines<CDim, Degree>::Impl`](https://ddc.mdls.fr/classddc_1_1UniformBSplines_1_1Impl.html) with arguments describing the uniform cells on which the splines are defined.

### Generating the Grid from the Spline Space

The grid is constructed from the Greville abcissae of the b-splines using DDC as follows:

```cpp
auto constexpr SplineXBoundary = ddc::BoundCond::PERIODIC;
auto constexpr SplineVxBoundary = ddc::BoundCond::HERMITE;

using SplineInterpPointsX
        = ddc::GrevilleInterpolationPoints<BSplinesX, SplineXBoundary, SplineXBoundary>;
using SplineInterpPointsVx
        = ddc::GrevilleInterpolationPoints<BSplinesVx, SplineVxBoundary, SplineVxBoundary>;

struct GridX : SplineInterpPointsX::interpolation_discrete_dimension_type
{
};
struct GridVx : SplineInterpPointsVx::interpolation_discrete_dimension_type
{
};
```

In the body of the simulation these objects must be initialised as global objects. This is done as follows:

```cpp
ddc::init_discrete_space<GridX>(SplineInterpPointsX::get_sampling<GridX>());
ddc::init_discrete_space<GridVx>(SplineInterpPointsVx::get_sampling<GridVx>());
```

The index range on which the simulation will be defined can then be constructed:

```cpp
IdxRange<X> idx_range_x(SplineInterpPointsX::get_domain<GridX>());
IdxRange<Vx> idx_range_vx(SplineInterpPointsVx::get_domain<GridVx>());
```

#### Spline Dependent Grids with Input

Spline-dependent grids are very common in Gyselalib++ but the required input for initialising these grids depends on:
- The uniformity of the cells on which the splines are defined
- The method of defining the grid points from the splines

This information is often provided as an input via PDI so Gyselalib++ also provides a helper function to initialise both the splines and the grids from a PDI configuration. This can be called as follows:

```cpp
IdxRange<X> const idx_range_x = init_spline_dependent_idx_range<
        GridX,
        BSplinesX,
        SplineInterpPointsX>(conf_gyselalibxx, "x");
IdxRange<Vx> const idx_range_vx = init_spline_dependent_idx_range<
        GridVx,
        BSplinesVx,
        SplineInterpPointsVx>(conf_gyselalibxx, "vx");
```

### Defining the Spline Interpolation

The spline interpolation requires a spline builder implementing the requested closure conditions, and a spline evaluator.

The spline builders are defined as:

```cpp
using SplineXBuilder = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX,
        GridX,
        SplineXBoundary,
        SplineXBoundary,
        ddc::SplineSolver::LAPACK>;

using SplineVxBuilder = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesVx,
        GridVx,
        SplineVxBoundary,
        SplineVxBoundary,
        ddc::SplineSolver::LAPACK>;
```

where `ddc::SplineSolver::LAPACK` indicates the solver that will be used to solve the matrix equation, and `Kokkos::DefaultExecutionSpace` indicates that the builder will act on the default execution space (GPU if enabled, for more information about execution spaces see the [Kokkos documentation](https://kokkos.org/kokkos-core-wiki/API/core/execution_spaces.html)).

In the body of the simulation an instance of the spline builder is created using the index range of the interpolation points:

```cpp
SplineXBuilder const builder_x(idx_range_x);
SplineVxBuilder const builder_vx(idx_range_vx);
```

The spline evaluators are defined as:

```cpp
using SplineXEvaluator = ddc::SplineEvaluator<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX,
        GridX,
        ddc::PeriodicExtrapolationRule<X>,
        ddc::PeriodicExtrapolationRule<X>>;

using SplineVxEvaluator = ddc::SplineEvaluator<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesVx,
        GridVx,
        ddc::ConstantExtrapolationRule<Vx>,
        ddc::ConstantExtrapolationRule<Vx>>
```

In the body of the simulation an instance of the spline evaluator operator is created from the extrapolation conditions:

```cpp
ddc::PeriodicExtrapolationRule<X> bv_x_min;
ddc::PeriodicExtrapolationRule<X> bv_x_max;

// The extrapolation is constant and has the same value as that of the spline
// at the first coordinate of the index range
ddc::ConstantExtrapolationRule<Vx> bv_vx_min(ddc::coordinate(idx_range_vx.front()));

// The extrapolation is constant and has the same value as that of the spline
// at the last coordinate of the index range
ddc::ConstantExtrapolationRule<Vx> bv_vx_max(ddc::coordinate(idx_range_vx.back()));

SplineXEvaluator const spline_x_evaluator(bv_x_min, bv_x_max);
SplineVxEvaluator const spline_vx_evaluator(bv_vx_min, bv_vx_max);
```

### Leveraging Common Geometries

In what we have presented so far we have described several structures which are specific to the geometry that we use in this simulation but are not specific to the simulation itself. In order to reduce code duplication and the time spent constructing such geometries, Gyselalib++ provides some files describing common geometries. These can be found in `src/geometry.../geometry/geometry.hpp`.

These files also contain several type aliases which simplify the notation of these types in the body of the code.

### Constructing the Semi-Lagrangian Advection Operators

Now we have defined both the splines and the grid we can construct our advection operators using spline interpolation:

```cpp
PreallocatableSplineInterpolator const spline_x_interpolator(builder_x, spline_x_evaluator);
PreallocatableSplineInterpolator const spline_vx_interpolator(builder_vx, spline_vx_evaluator);

BslAdvectionSpatial<GeometryXVx, GridX> const advection_x(spline_x_interpolator);
BslAdvectionVelocity<GeometryXVx, GridVx> const advection_vx(spline_vx_interpolator);
```

### Implementing the Operator Splitting Strategy

The overall advection operator is constructed using an operator implementing Strang splitting. For the XVx geometry Gyselalib++ provides the `SplitVlasovSolver` operator for this:

```cpp
SplitVlasovSolver const vlasov(advection_x, advection_vx);
```

### Solving the Poisson Equation

The Poisson equation is solved in two stages:
1. Evaluate the charge density
2. Solve the equation using an FFT solver

#### Evaluating the Charge Density

The quadrature scheme for the charge density operator first requires the definition of the quadrature coefficients using the chosen scheme:

```cpp
DFieldMemVx const quadrature_coeffs(neumann_spline_quadrature_coefficients<
                                    Kokkos::DefaultExecutionSpace>(mesh_vx, builder_vx));
```

Here `DFieldMemVx` is a field of doubles defined along the Vx dimension. The keyword `Mem` indicates that this object allocates memory on instantiation. It should not be used directly, the field must always be extracted. For more explanations see [Using DDC in Gyselalib++](./DDC_in_gyselalibxx.md).

The charge density calculator can be expressed simply using the `Quadrature` class. Gyselalib++ provides an operator for the XVx case which can be used as an example for other geometries.
This operator is assembled from this field of quadrature coefficients as follows:

```cpp
ChargeDensityCalculator rhs(get_field(quadrature_coeffs));
```

#### Solving the Equation with an FFT Solver

The FFT solver is constructed using the index range on which it operates:

```cpp
FFTPoissonSolver<IdxRangeX> fft_poisson_solver(idx_range_x);
```

This operator and the charge density calculator are bound together using a class (specific to the XVx geometry) called `QNSolver`:

```cpp
QNSolver const poisson(fft_poisson_solver, rhs);
```

### Assembling the Time-Stepping Scheme

Finally the predictor-corrector (specific to the XVx geometry) is constructed from the Vlasov operator and the Poisson operator:

```cpp
PredCorr const predcorr(vlasov, poisson);
```

---

## Initial Conditions
