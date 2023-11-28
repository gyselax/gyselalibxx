# Poisson Solver

The Poisson solver is designed to solve the following Poisson equation:

$$ -\Delta \phi(x) = \sum_s \int q_s f_s(x,y,v_x,v_y) dv_x dv_y $$

This calculation is split into two parts. Firstly the charge density is calculated:

$$ \rho_{q,s}(x, y) = \int q_s f_s(x,y,v_x,v_y) dv_x v_y $$

Then the basic Poisson equation:

$$ -\Delta \phi(x) = \sum_s \rho_{q,s}(x,y) $$

is solved.

## Charge Density

The charge density is calculated by integrating the distribution function. Currently this is done by creating an intermediate representation of a spline and integrating that representation. This is costly and should be replaced by a ChargeDensityCalculator which uses a spline quadrature.

## Poisson Solver

The Poisson equation can be solved with a variety of different methods. Here we have implemented:

-   FftPoissonSolver

These classes return the electric potential $\phi$ and the electric field $\frac{d \phi}{dx}$.

The FftPoissonSolver does not calculate the electric field using the Fourier modes. Rather it uses a spline interpolation to approximate this value. This interpolation is calculated by the operator ElectricField.

