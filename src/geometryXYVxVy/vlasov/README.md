# Vlasov solver

This folder contains methods for solving the Vlasov equation:
```math
\partial_t f(t, x, v) + v_x \partial_x f + v_y \partial_y f + \frac{E_x}{m} \partial_{v_x} f + \frac{E_y}{m} \partial_{v_x} f = 0
```

Where $f$ is the distribution function, $x$ and $y$ are the spatial variables, \(v_x\) and \(v_y\) are the velocity variables, \(E_x\) and \(E_y\) are the components of the electric field and $m$ is the mass of the considered plasma species.

The implemented solvers are:
- SplitVlasovSolver : Solves the Vlasov equation using Strang splitting
- MpiSplitVlasovSolver : Solves the Vlasov equation using Strang splitting and MPI transposes between a X2Dsplit and a V2Dsplit layout.
