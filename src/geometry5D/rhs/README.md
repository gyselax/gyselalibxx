# RHS 

The `rhs` folder contains any operator that appears on the right-hand-side (rhs) of a Boltzmann equation. These operators are often referred to as `sources`, to differentiate them from the advective terms that appear on the left-hand-side of a Boltzmann equation. The term `sources` does not necessarily imply that the considered operator injects particles in the plasma. It can act as a source or a sink of particles, or as a collision term that conserves the number of particles for instance.

The currently implemented sources are:
- Collisions : the collision (vpar,mu) operator (translate from Fortran version in C++ in koliop submodule)
