# Boltzmann solver

The `boltzmann` folder contains methods for solving a Boltzmann equation. Such methods typically take the distribution function and the electric field computed at a given time $t$, and return the value of the distribution function at a time $t+dt$, where $dt$ is the timestep of the simulation. A Boltzmann equation refers to an advection equation in phase space with sources. In the simplified 1D geometry in space and velocity it has the general form

$`\partial_t f + v \partial_x f + \frac{E}{m}\, \partial_v f = S(f)`$

Where $f$ is the distribution function, $x$ and $v$ are the space and velocity variables respectively, $E$ is the electric field and $m$ is the mass of the considered plasma species.

The operator \( S(f) \) represents any source terms.
Different types of source terms include:

- [Kinetic sources](../../../docs/latex/geometryXVx/rhs/kinetic_source.pdf)
- [Krook sources](../../../docs/latex/geometryXVx/rhs/krook_source.pdf)

To incorporate spatial dependencies in source terms, mask functions are often used.
These functions take the value 1 where the source is active and 0 otherwise.
An example of such a mask function can be found in
[Mask with tanh](../../../docs/latex/geometryXVx/rhs/mask_tanh.pdf).

When setting \( S(f) = 0 \), the equation is commonly referred to as the Vlasov equation.

The implemented Boltzmann solvers are:

- SplitRightHandSideSolver
- SplitVlasovSolver
