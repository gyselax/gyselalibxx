# PDE Solvers

This folder contains all generalised code describing different methods for solving PDEs (Partial Differential Equations). Solvers exist for the following PDEs:

- Poisson's equation : $-\Delta \phi = \rho$

## Poisson's equation

The following methods exist for solving Poisson's equations:
- FFTPoissonSolver

These solvers implement 2 interfaces:
- `chunk_span_type operator()(chunk_span_type phi, chunk_span_type rho) const`
- `chunk_span_type operator()(chunk_span_type phi, vector_span_type E, chunk_span_type rho) const`

The second interface calculates $\phi$ the solution to the equation but also $E = - \nabla \phi$.

