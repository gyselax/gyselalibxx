# Sheath simulations

## Description 
The `sheath` executable allows the study of kinetic plasma-wall interaction. A Boltzmann-Poisson system is solved for the electric field and the distribution function for both electrons and ions species. The model is one dimensional in space and velocity (1D1V).   

The particularity of the model lies in the description of the wall region, which is treated using immersed boundary condition (or penalization). Electrons and ions are absorbed inside this region. The inertia difference between the two species drives the appearance of a positively charged layer (the sheath) at the plasma boundary. More details about the sheath physics along with the need for non-equidistant meshes to conduct such simulations can be found in [1,2].  

## Usage
After building the code, run the executable located in `build/simulations/sheath/`. To use the default simulation parameters the user can provide the `--dump-config` option when launching the executable.

## Recommended parameters
Relevant parameters for a test-case simulation are given in the `sheath.yaml.hpp` file. To observe a steady state plasma-wall interaction with interesting physical features (supersonic flow of ions, truncation of fast velocity electrons, etc.) see the paramerers of [2].

## References
- [1] E. Bourne, Y. Munschy, V. Grandgirard, M. Mehrenberger, and P. Ghendrih, Non-Uniform Splines for Semi-Lagrangian Kinetic Simulations of the Plasma Sheath (2022)
- [2] Y. Munschy, E. Bourne, P. Ghendrih, G. Dif-Pradalier, Y. Sarazin, V. Grandgirard, and P. Donnel, Kinetic plasma-wall interaction using immersed boundary conditions (2023)
