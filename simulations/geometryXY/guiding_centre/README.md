# Guiding centre (X,Y) simulation

## Equations

The `guiding_centre` executable deals with the following model:

```math
    \partial_t f(t, x, y) + (E(t, x, y)\wedge e_z) \cdot \nabla f(t, x, y) = 0, \\
    - \Delta \phi(t, x, y)  = -\partial_x^2 \phi(t, x, y) -\partial_y^2 \phi(t, x, y) = f(t, x, y), \\
    E(t, x, y) =  - \nabla \phi(t, x, y).
```

**Remark:** the advection field of the advection equation is given by $`\frac{E\wedge B}{|B|^2} = E(t, x, y)\wedge e_z`$ with $`B`$ the magnetic field supposed of norm 1, (and $`e_z`$ a vector perpendicular to the $`(x,y)`$ plane in the direct direction).

### Operators

The simulations uses the following operators:

- advection equation: BslAdvection1D operator with a Strang splitting along $`x`$ and $`y`$.
The time integration methods applied to solve the characteristic equation are explicit Euler methods;
- Poisson equation: FFTPoissonSolver solver using FFT to solve the Poisson equation on a periodic domain (and compute the electric field);
- equations coupling: PredCorrRK2XY using a RK2 time integration method.

### Execution space

The simulation runs on GPU/CPU.

### Output files

The `guiding_centre` executable creates in the working folder an `output/` folder containing the output files: density, electrostatic potential and electric field values on the grid for some regular time steps.

### Command line

```shell
./guiding_centre_XY <path_to_params.yaml>
```

## Simulation

### Grid

The simulation runs on a grid on $`[0, 4\pi]\times[0, 2\pi]`$  with 64 uniformly distributed points in each direction.

### Test case: Kelvin-Helmholtz instability test case

The chosen initial conditions are

```math
    f(0, x, y) = f_{\text{eq}}(x,y) + \varepsilon\cos(kx),  \\
    f_{\text{eq}}(x,y) = \sin(y)
```

with $`\varepsilon = 0.015`$ the amplitude of perturbation and $`k = 2\pi/ L_x = 0.5`$ the mode of the perturbation.

See more in [initialisations](./../../../src/geometryXY/initialisation/README.md).

### Results treatments - diagnosis

In `post-process/PythonScripts/geometryXY/` folder, python files are available to plot the results saved in the `output/` folder:

- `animation.py` displays and can save the curves of the density $`f`$ and the electrostatic potential $`\phi`$ for each time steps saved in the `output/` folder.

The command line is given by

```shell
python3 animation.py --name=<name_file_to_save> --folder=<path_to_output>
```

- `energy_conservation.py` plots the evolution of the energy of the system along the simulation:

```math
d\mathcal{M}: t\mapsto \int_{\Omega} |E(t)|^2 - \int_{\Omega} |E_0|^2
```

The command line is given by

```shell
python3 energy_conservation.py --name=<name_file_to_save> --folder=<path_to_output>
```

- `mass_conservation.py` plots the evolution of the mass of the system along the simulation:

```math
d\mathcal{W}: t\mapsto \int_{\Omega} f(t) - \int_{\Omega} f_0
```

The command line is given by

```shell
python3 mass_conservation.py --name=<name_file_to_save> --folder=<path_to_output>
```

- `plot_individual_curve.py` plots the curves of the density $`f`$ and the electrostatic potential $`\phi`$ for a given time

The command line is given by

```shell
python3 plot_individual_curve.py --name=<name_file_to_save> --time=<given_time> --folder=<path_to_output>
```

- `plot_L2_norms.py` plots the $`\mathcal{L}^2`$ norms of the density $`f - f_{eq}`$ and the electrostatic potential $`\phi`$ along the simulation

The command line is given by

```shell
python3 plot_L2_norms.py --name=<name_file_to_save> --folder=<path_to_output>
```

## Contents

- `guiding_centre.cpp`: executable of a guiding-centre equation on $`(x,y)`$ geometry with Kelvin-Helmholtz instability test case initial conditions.
- `params.yaml`: contains the parameters of the simulation. It needs to be added at the command line to launch the executable.

## References

[1]     Eric Sonnendrücker, Jean Roche, Pierre Bertrand and Alain Ghizzo. "The Semi-Lagrangian Method for the Numerical Resolution of the Vlasov Equation", on Journal of Computational Physics, Oct. 1998.

[2]     Nicolas Crouseilles, Michel Mehrenberger, Eric Sonnendrücker. "Conservative semi-Lagrangian schemes for Vlasov equations". Oct. 2009.
