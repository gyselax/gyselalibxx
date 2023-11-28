# Poisson Solver

The Poisson solver is designed to solve the following Poisson equation:

$$ -\frac{d^2 \phi}{d x^2}(x) = \sum_s \int_v q_s f_s(x,v) dv $$

This calculation is split into two parts. Firstly the charge density is calculated:

$$ \rho_{q,s}(x) = \int_v q_s f_s(x,v) dv $$

Then the basic Poisson equation:

$$ -\frac{d^2 \phi}{d x^2}(x) = \sum_s \rho_{q,s}(x) $$

is solved.

## Charge Density

The charge density is calculated by integrating the distribution function. The simplest way of doing this is using the ChargeDensityCalculator class which takes a quadrature method.

For historical reasons the class SplineChargeDensityCalculator also exists. This class creates an intermediate representation of a spline and integrates that representation. This is equivalent to using ChargeDensityCalculator and providing a spline quadrature. The quadrature method should be preferred as the numerical cost is significantly smaller.

## Poisson Solver

The Poisson equation can be solved with a variety of different methods. Here we have implemented:

-   FftPoissonSolver
-   FemNonPeriodicPoissonSolver
-   FemPeriodicPoissonSolver

These classes return the electric potential $\phi$ and the electric field $\frac{d \phi}{dx}$.

The FftPoissonSolver does not calculate the electric field using the Fourier modes. Rather it uses a spline interpolation to approximate this value. This interpolation is calculated by the operator ElectricField.
