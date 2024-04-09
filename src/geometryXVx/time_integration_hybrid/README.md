# Time integration hybrid

Code for time integrators that solves a Boltzmann-Poisson system for the electrostatic potential and the electron and ion distribution function, along with a fluid model to describe a fluid species. A predictor-corrector scheme is used to solve the Vlasov Poisson system coupled to the fluid equation on one timestep $\Delta t$ as
1. Find the electrostatic potential at time $t+\Delta t/2$ following: 
- Solve Boltzmann equation for the distribution function at time $t+\Delta t/2$; 
- Using this distribution function, solve Poisson equation for the electrostatic potential at time $t+\Delta t/2$.
2. With the potential computed at time $t+\Delta t/2$:
- Solve Boltzmann equation for the distribution function at time $t+\Delta t$; 
- Using the distribution function at time $t+\Delta t$, solve Poisson equation for the electrostatic potential at time $t+\Delta t$;
- Solve the fluid model for the fluid species fluid moments at time $t+\Delta t$.

For more information about the Boltzmann-Poisson time integrator, see [time\_integration](./../time_integration/README.md).  

The implemented time integrators that take into account fluid species are: 
- PredCorrHybrid