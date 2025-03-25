# Geometry (x, v\_x)

The `geometryXVx` folder contains all the code describing methods which are specific to a geometry with 1 spatial dimension and 1 velocity dimension. It is broken up into the following sub-folders:

- [boltzmann](./boltzmann/README.md) : Solvers for a Boltzmann equation.
- [geometry](./geometry/README.md) : All the dimension tags used for a simulation in the geometry.
- [initialisation](./initialisation/README.md) : Initialisation methods for the distribution function.
- [poisson](./poisson/README.md) : Code describing the Quasi-Neutrality solver.
- [rhs](./rhs/README.md) : Code describing the operators on the right hand side of the Boltzmann equation; namely sources, sinks and collisions.
- [time\_integration](./time_integration/README.md) : Time integrators for a Boltzmann-Poisson system of equations.
- [utils](./utils/README.md) : Miscellaneous utility functions.
