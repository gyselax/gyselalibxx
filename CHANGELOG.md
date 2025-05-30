# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [UNRELEASED]

### Added

### Fixed

- Fix paths in root `CMakeLists.txt` file to ensure it can be correctly used in a submodule.
- Update remaining use of `ddc::Coordinate` to use Gyselalib++ conventions (`Coord`).
- Update coding conventions to match what is applied.

### Changed

### Deprecated

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
