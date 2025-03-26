# Initialisation on (x,y) geometry

Describes different initial conditions for the simulations in $`(x,y)`$ geometry.

(See more on the simulation in [simulations](./../../../simulations/geometryXY/README.md).)

## Kelvin-Helmholtz instability test case

For this test case, the KelvinHelmholtzInstabilityInitialisation sets the initial conditions at

```math
    f(0, x, y) = f_{\text{eq}}(x,y) + \varepsilon\cos(kx),  \\
    f_{\text{eq}}(x,y) = \sin(y)
```

with $`\varepsilon = 0.015`$ the amplitude of perturbation and $`k = 2\pi/ L_x = 0.5`$ the mode of the perturbation.

We suppose the domain periodic on $`x`$ and $`y`$.

These are applied in the guiding-centre equations. See in [guiding-centre](./../../../simulations/geometryXY/guiding_centre/README.md).
