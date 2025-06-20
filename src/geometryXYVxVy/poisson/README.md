# Quasi-Neutrality Solver

The Quasi-Neutrality solver is designed to solve the following Quasi-Neutrality equation:

$$ -\Delta \phi(x) = \sum_s \int q_s f_s(x,y,v_x,v_y) dv_x dv_y $$

This calculation is split into two parts. Firstly the charge density is calculated:

$$ \rho_{q,s}(x, y) = \int q_s f_s(x,y,v_x,v_y) dv_x v_y $$

Then the basic Poisson equation:

$$ -\Delta \phi(x) = \sum_s \rho_{q,s}(x,y) $$

is solved.

## Charge Density

The charge density is calculated by integrating the distribution function.

## Quasi-Neutrality Solver

The Quasi-Neutrality equation can be solved with a variety of different methods by combining Poisson solvers and charge density solvers.

These classes return the electric potential $\phi$ and the electric field $\frac{d \phi}{dx}$.
