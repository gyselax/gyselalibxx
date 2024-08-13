# PDE Solvers

This folder contains all generalised code describing different methods for solving PDEs (Partial Differential Equations). Solvers exist for the following PDEs:

- Poisson's equation : $-\Delta \phi = \rho$

## Poisson's equation

The following methods exist for solving Poisson's equations:
- FFTPoissonSolver

These solvers implement 2 interfaces:
- `field_type operator()(field_type phi, field_type rho) const`
- `field_type operator()(field_type phi, vector_field_type E, chunk_field_type rho) const`

The second interface calculates $\phi$ the solution to the equation but also $E = - \nabla \phi$.

