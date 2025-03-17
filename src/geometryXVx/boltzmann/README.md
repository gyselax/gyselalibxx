# Boltzmann solver

The `boltzmann` folder contains methods for solving a Boltzmann equation. Such methods typically take the distribution function and the electric field computed at a given time $t$, and return the value of the distribution function at a time $t+dt$, where $dt$ is the timestep of the simulation. A Boltzmann equation refers to an advection equation in phase space with sources. In the simplified 1D geometry in space and velocity it has the general form 

\(\partial_t f + v \partial_x f + \frac{E}{m}\, \partial_v f = S(f)\)

Where $f$ is the distribution function, $x$ and $v$ are the space and velocity variables respectively, $E$ is the electric field and $m$ is the mass of the considered plasma species. The $S(f)$ operator refers to any source terms (including collisions). When taking $S(f)=0$ this equation is often referred to as a Vlasov equation.

The implemented Boltzmann solvers are: 
- SplitRightHandSideSolver
- SplitVlasovSolver