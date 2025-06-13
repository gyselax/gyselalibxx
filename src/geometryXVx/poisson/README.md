# Quasi-Neutrality Solver

The Quasi-Neutrality solver is designed to solve the following Quasi-Neutrality equation:

$$ -\frac{d^2 \phi}{d x^2}(x) = \sum_s \int_v q_s f_s(x,v) dv $$

This calculation is split into two parts. Firstly the charge density is calculated:

$$ \rho_{s}(x) = \sum_s \int_v q_s f_s(x,v) dv $$

Then the basic Poisson equation:

$$ -\frac{d^2 \phi}{d x^2}(x) = \rho_{q,s}(x) $$

is solved.

## Charge Density

The charge density is calculated by integrating the distribution function. The simplest way of doing this is using the ChargeDensityCalculator class which takes a quadrature method.

## Poisson Solver

The Quasi-Neutrality equation can be solved with a variety of different methods. Here we have implemented:

- FftQNSolver
- FemNonPeriodicQNSolver
- FemPeriodicQNSolver

These classes return the electric potential $\phi$ and the electric field $\frac{d \phi}{dx}$.
