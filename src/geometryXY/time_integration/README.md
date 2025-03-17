# Predictor-corrector methods

## Predictor-corrector based on RK2

The PredCorrRK2XY defines a predictor-corrector based on the RK2 time integration method. 
The RK2 method is implemented in the ITimeStepper method (see [rk2](../../timestepper/README.md)). 

This time integration method is applied to a a guiding-centre equations system documented in [guiding\_centre](./../../../simulations/geometryXY/guiding_centre/README.md). 
```math
    \partial_t f(t, x, y) + (E(t, x, y)\wedge e_z) \cdot \nabla f(t, x, y) = 0, \\
    - \Delta \phi(t, x, y)  = -\partial_x^2 \phi(t, x, y) -\partial_y^2 \phi(t, x, y) = f(t, x, y), \\
    E(t, x, y) =  - \nabla \phi(t, x, y).
```


Applied to this equations system, we have, writing \(\{t^n\}_n \text{ the time discretisation, and } f^n = f(t^n, \cdot, \cdot)\):


For $n\geq 0$,

* Advect on a half time step:

 1. From $f^n$, we compute $\phi^n$ with the Poisson-like equation (FFTPoissonSolver);

 2. From $\phi^n$, we compute $E^n$ by deriving (FFTPoissonSolver);

 3. From $f^n \text{ and } E^n$, we compute $f^{n+1/2}$ by advecting (BslAdvection1D) on $\frac{dt}{2}$;

* Advect on a full time step:

 4. From $f^{n+1/2}$, we compute $\phi^{n+1/2}$ with the Poisson-like equation (FFTPoissonSolver);

 5. From $\phi^{n+1/2}$, we compute $E^{n+1/2}$ by deriving (FFTPoissonSolver);

 6. From $f^n \text{ and } E^{n+1/2}$, we compute $f^{n+1}$ by advecting (BslAdvection1D) on $dt$.

