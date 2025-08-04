# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [UNRELEASED]

### Added

- Port `PolarSplineEvaluator` methods to GPU.
- Add methods to `PolarSplineEvaluator` to avoid unnecessary creation of fields of coordinates.

### Fixed

- Modify `ruche.v100/environment.sh` file to fix tokamaxi simulation segfault issues.
- Fix type of derivatives stored in `DerivFieldMem` and `DerivField` types.
- Fixed memory error in `PolarSplineFEMPoissonLikeSolver`.
- Remove non-parallelisable loop in `PolarSplineFEMPoissonLikeSolver::init_nnz_per_line`.
- Remove use of `std::cyl_bessel_j` which is not available in libc++.
- Fix `mi250.hipcc.adastra.spack` toolchain.
- Fix uninitialised values being used as an initial guess for the result of the matrix equation in `PolarSplineFEMPoissonLikeSolver`.
- Fix missing grids when calling `collect_grids_on_dim_t`.

### Changed

- Change interface of `EdgeTransformation::search_for_match` to return an `out_of_bounds_idx` instead of a boolean.
- Change spack setup in CPU installation script (`prepare.sh`) to create and use independent spack installation.
- Change template parameters of `PolarSplineEvaluator` to add execution and memory space information.
- Allow `get_idx_range` to be called from a GPU execution space.
- Uniformise toolchains.
- Allow batch CSR convergence parameters to be specified in the constructor of `PolarSplineFEMPoissonLikeSolver`.
- Change the internals of `PolarSplineFEMPoissonLikeSolver` to precalculate fewer values.
- Change the internals of `PolarSplineFEMPoissonLikeSolver` to avoid calls to DDC's internals.
- Clean up code in `BslImplicitPredCorrRTheta::operator()`.

### Deprecated

### Removed

- Remove deprecated method `PolarBSplines::integrals`.
- Remove unhelpful `PolarSpline` classes in favour of `DField<IdxRange<PolarBSplines>>` types.
- Remove unused broken method `PolarSplineEvaluator::integrate`.

## [v0.2.0] - 2025-07-03

### Added

- Add a Gyroaverage operator with tests for circular geometry.
- Curvilinear coordinate change classes have an O-point method to retrieve the O-point in the non-curvilinear coordinates.
- Add a batched `operator()` to `DiscreteToCartesian` allowing a field of coordinates to be converted.
- Add a `LiePoissonBracket::operator()` overload which takes a 2D tensor as the second argument to the bracket.
- Add a function `PDI_expose_vector_field` to output a vector field via PDI.
- Add a `control_points` method to `DiscreteToCartesian` to allow all control points to be retrieved at once.

### Fixed

- Fix uninitialized warning in the `Tensor` class.
- Fix unused `m_magnetic_field` variable in `MaxwellianEquilibrium` class.
- Fix break points incorrectly labelled as knots.
- Fix minimum version requirement of Kokkos.
- Fix tolerance of floating point comparisons in JacobianMatrixAndJacobianCoefficients and MultipatchSplineEvaluatorTest tests
- Fix unnecessary `std::move` calls.
- Fix missing assertion in `LeviCivitaTensor` to prevent division by 0 when Jacobian is calculated at singular point.
- Fix bad result of `is_tensor_v` for `IdentityTensor`.

### Changed

- Bumped the minimum CMake version to 3.25 to benefit from the `add_subdirectory(.. SYSTEM)` feature.
- Change the interface of `IVlasovSolver` and `IQNSolver` in `geometryXYVxVy` to store the electric field in a `VectorField`.
- Integration of `ddc::StridedDiscreteDomain` by making `IdxRangeSlice` a type alias.
- The parameter `iter_start` has been removed from the constructor of `RestartInitialisation`.
- Generalise `compute_coeffs_on_mapping` to work with any mapping.
- Rely on GPU-aware MPI to allow GPU-direct MPI for `MPITransposeAllToAll`.
- Curvilinear coordinate change classes take a `Coord` type to specify the O-point in the constructor.
- Allow `init_discrete_space` to be used to initialise `PolarBSplines` with a GPU-based `DiscreteToCartesian` coordinate change operator.

### Deprecated

## [v0.1.1] - 2025-05-30

### Fixed

- Fix paths in root `CMakeLists.txt` file to ensure it can be correctly used in a submodule.
- Update remaining use of `ddc::Coordinate` to use Gyselalib++ conventions (`Coord`).
- Update coding conventions to match what is applied.

## [v0.1.0] - 2025-05-28

### Added

- First release of Gyselalib++
- Advection operators
  - 1D Semi-Lagrangian spatial advection ($` \frac{df_s}{dt}= \sqrt{\frac{m_e}{m_s}} v \frac{\partial f_s}{\partial x} `$)
  - 1D Semi-Lagrangian velocity advection ($` \frac{df_s}{dt}= q_s \sqrt{\frac{m_e}{m_s}} E \frac{\partial f_s}{\partial v} `$)
  - 1D Semi-Lagrangian advection with a provided advection field
  - 2D Semi-Lagrangian advection on a polar plane with a provided advection field
- Collisions
  - Collision operator in $`(v_\parallel,\mu)`$
- Coordinate transformation operators and tools
  - Coordinate transformation operators
    - Triangular Barycentric coordinates <-> Cartesian coordinates
    - Circular coordinates <-> Cartesian coordinates
    - Cylindrical coordinates <-> Cartesian coordinates
    - Tokamak-shaped Czarny coordinates <-> Cartesian coordinates
    - Toroidal coordinates -> Cylindrical coordinates
    - Discrete coordinates -> Cartesian coordinates
    - Identity transformation
    - Composite coordinate transformation
  - Tools to manage coordinate transformations by:
    - Getting the inverse Jacobian matrix at a given coordinate
    - Getting the inverse Jacobian matrix at the O-point (to provide explicit equations without an if)
    - Evaluate the metric tensor at a given coordinate
    - Map a vector from one vector space to another
- Additional data types
  - DerivativeField to store a field and its boundary derivatives
  - VectorField
  - Tensor type and tools
    - Levi-Civita tensor
    - Identity tensor
    - Tensor multiplication operator
- Interpolation operators
  - Lagrange interpolation
  - Spline interpolation
  - Polar spline evaluation
- General Mathematical tools
  - Methods for calculating the L-norms
  - Derivative calculators
    - Finite differences method (with and without known boundary values)
    - Derivatives from 1D or 2D spline representations
    - Constant derivatives of a known value
  - Miscellaneous
    - sum
    - norm
    - modulo
    - pow for integer powers
    - factorial
    - min
    - max
    - determinant
    - inverse
    - scalar product
    - tensor product
- Solvers for matrix equations with different sparsity patterns
  - Banded matrix
  - Batched CSR (compressed sparse row) matrix
  - Batched ELL matrix
  - Batched tridiagonal matrix
  - Matrix with dense corners
  - Matrix with dense bottom right-hand corner
  - Dense matrix
  - Positive-definite symmetric tridiagonal matrix
  - Periodic banded matrix (banded plus top-right and bottom-left corners)
- MPI parallelisation tools
  - MPI all to all transposition
- Solvers for partial differential equation (PDEs)
  - 1D Finite elements method (FEM)
  - 1D Fast Fourier transform (FFT)
  - Polar Poisson-like solver to solve $` - \nabla \cdot (\alpha \nabla \phi) + \beta \phi = \rho `$ on a polar domain
- Quadrature methods and tools
  - Definition of quadrature coefficients
    - Gauss-Legendre
    - Trapezoid
    - Simpson
    - Spline-based quadrature (with homogeneous Neumann boundary conditions or with an equal number of spline bases and interpolation points)
  - ND coefficients from multiple 1D quadrature coefficients
  - Definition of quadrature coefficients to calculate a volume on general coordinates
- Timestepping methods
  - Euler
  - Crank-Nicolson
  - 2nd order Runge-Kutta (RK2)
  - 3rd order Runge-Kutta (RK3)
  - 4th order Runge-Kutta (RK4)

### Notes

- This is an early development release (`v0.x`); APIs may change
