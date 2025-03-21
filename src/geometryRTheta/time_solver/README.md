# Predictor-corrector methods

This folder contains different predictor-corrector methods:

- Predictor-corrector: BslPredCorrRP;
- Second order explicit predictor-corrector: BslExplicitPredCorrRP;
- Second order implicit predictor-corrector: BslImplicitPredCorrRP;

The second order explicit predictor-corrector and the second order implicit predictor-corrector are detailed in Edoardo Zoni's article [1].

The studied equations system is Vlasov-Poisson equations

```math
 L\phi = - \nabla \cdot (\alpha \nabla \phi) + \beta \phi = \rho,
\\ E = - \nabla \phi,
\\ (A_x, A_y) = (-E_y, E_x),
\\ \partial_t \rho + A_x \partial_x \rho +  A_y \partial_y \rho = 0,
\\ \rho(t = 0, r, \theta) = \rho_0(r,\theta).
```

Let's write $`\{t^n\}_n \text{ the time discretisation, and } f^n = f(t^n, \cdot, \cdot)`$.

## Predictor-corrector

This method is a Runge-Kutta 2 method: for $n\geq 0$,

Advect on a half time step:

 1. From $\rho^n$, we compute $\phi^n$ with the Poisson-like equation (IQNSolver);

 2. From $\phi^n$, we compute $E^n$ by deriving (IQNSolver);

 3. From $\rho^n \text{ and } A^n$, we compute $\rho^{n+1/2}$ by advecting (IAdvectionRTheta) on $\frac{dt}{2}$ with one of the selected time integration methods (ITimeStepper);

Advect on a full time step:

 4. From $\rho^{n+1/2}$, we compute $\phi^{n+1/2}$ with the Poisson-like equation (IQNSolver);

 5. From $\phi^{n+1/2}$, we compute $E^{n+1/2}$ by deriving (IQNSolver);

 6. From $\rho^n \text{ and } A^{n+1/2}$, we compute $\rho^{n+1}$ by advecting (IAdvectionRTheta) on $dt$ with one of the selected time integration methods (ITimeStepper).

## Explicit predictor-corrector

This method is a second order explicit predictor-corrector: for $n\geq 0$,

Prediction:

 1. From $\rho^n$, we compute $\phi^n$ with the Poisson-like equation (IQNSolver);

 2. From $\phi^n$, we compute $E^n$ by deriving (IQNSolver);

 3. From $\rho^n \text{ and } A^n$, we compute $\rho^P$ by advecting (IAdvectionRTheta) on $dt$;
  - We write $X^P$ the characteristic feet such that $`\partial_t X^P = A^n(X^n)`$.

Correction:

 4. From $\rho^{P}$, we compute $\phi^{P}$ with the Poisson-like equation (IQNSolver);

 5. From $\phi^{P}$, we compute $E^{P}$ by deriving (IQNSolver);

 6. From $\rho^n \text{ and } \frac{A^{P}(X^n) + A^n(X^P)}{2}$, we compute $\rho^{C} = \rho^{n+1}$ by advecting (IAdvectionRTheta) on $dt$.

## Implicit predictor-corrector

This method is a second order implicit predictor-corrector: for $n\geq 0$,

Prediction:

 1. From $\rho^n$, we compute $\phi^n$ with the Poisson-like equation (IQNSolver);

 2. From $\phi^n$, we compute $E^n$ by deriving (IQNSolver);

 3. From $\rho^n \text{ and } A^n$, we compute $\rho^P$ by advecting (IAdvectionRTheta) on $\frac{dt}{2}$;
  - We write $X^P$ the characteristic feet such that it is the result of the implicit method:

```math
\partial_t X^k = \frac{A^n(X^n) + A^n(X^{k-1})}{2},  \qquad  X^k = X^n - \frac{dt}{2} \partial_t X^k.
```

Correction:

 4. From $\rho^{P}$, we compute $\phi^{P}$ with the Poisson-like equation (IQNSolver);

 5. From $\phi^{P}$, we compute $E^{P}$ by deriving (IQNSolver);

 6. From $\rho^n \text{ and } A^{P}$, we compute $\rho^{C} = \rho^{n+1}$ by advecting (IAdvectionRTheta) on $dt$.
   - We write $X^C$ the characteristic feet such that it is the result of the implicit method:

```math
\partial_t X^k = \frac{A^P(X^n) + A^P(X^{k-1})}{2},  \qquad  X^k = X^n - dt \partial_t X^k.
```

## References

[1] Edoardo Zoni, Yaman Güçlü, "Solving hyperbolic-elliptic problems on singular mapped disk-like domains with the
method of characteristics and spline finite elements", <https://doi.org/10.1016/j.jcp.2019.108889>, Journal of Computational Physics, 2019.

## Contents

- itimesolver.hpp: base class for the time solvers: ITimeSolverRP;
- bsl\_predcorr\_second\_order\_explicit.hpp: the second order explicit predictor-corrector (BslExplicitPredCorrRP);
- bsl\_predcorr\_second\_order\_implicit.hpp: the second order implicit predictor-corrector (BslImplicitPredCorrRP);
- bsl\_predcorr.hpp: the Runge-Kutta 2 method (BslPredCorrRP).
