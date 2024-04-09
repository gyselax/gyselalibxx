# Fluid solvers (MX)

The `fluidsolver` folder contains code that allows solving for models that describe a species considered as fluid. Such a fluid species is described by moments of the distribution function (density, particle flux, energy...) rather than by a complete distribution function. Equations describing the time evolution of such a fluid moments are called *fluid equations*. Such fluid equations describe the conservation of fluid moments: conservation of density, particle flux, energy, etc.

The currently implemented solvers are : 
- NullFluidSolver