# Neutral simulations

## Description 
The `neutral` executable allows the study of plasma-neutral interactions on a single magnetic field line. The magnetic field line is connected at both of its ends to a solid wall. A Boltzmann-Poisson system is solved for the electric field and the distribution function for both electrons and ions species. The model is one dimensional in space and velocity (1D1V).   

The wall region is treated using immersed boundary condition (or penalization). Electrons and ions are absorbed inside this region. The inertia difference between the two species drives the appearance of a positively charged layer (the sheath) at the plasma boundary. The neutrals are treated using a simple fluid model (advective, pressure-diffusion, etc.)

## Usage
After building the code, run the executable located in `build/simulations/neutral/`. To use the default simulation parameters the user can provide the `--dump-config` option when launching the executable.