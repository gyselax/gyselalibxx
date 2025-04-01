# Time Stepping Methods

Time stepping methods are methods for calculating how a value (or dimensioned values) evolves over time. Such an evolution should be expressed as a system of equations of the form:

```math
d\overrightarrow{x}/dt = f(\overrightarrow{x})\\
\overrightarrow{x}(t_0) = \overrightarrow{x}_0
```

The implemented methods are:

- Explicit Euler (first order) (Euler)
- Crank-Nicolson (second order) (CrankNicolson)
- Second order Runge Kutta (RK2)
- Third order Runge Kutta (RK3)
- Fourth order Runge Kutta (RK4)

These classes all contain an `update` method which carries out one time step of the algorithm.

## Explicit Euler method

- Scheme:
$x^{n+1} = x^{n} + dt f(x^{n+1})$

- Convergence order : 1.

## Crank-Nicolson method

- Scheme:
$x^{k} = x^{n} + \frac{dt}{2} \left( f(x^{n}) + f(x^k) \right)$ and
$x^{n+1} = x^{k+1}$ once converged.

- Convergence order : 2.

## RK2 method

- Scheme:
$`x^{n+1} = x^{n} + dt k_2`$
  - with
    - $`k_1 =  f(x^{n})`$,
    - $`k_2 =  f(x^{n} + \frac{dt}{2} k_1)`$,

- Convergence order : 3.

## RK3 method

- Scheme:
$`x^n = x^{n+1} - \frac{dt}{6}  \left( k_1 + 4 k_2 + k_3 \right)`$
  - with
    - $`k_1 =  f(x^{n})`$,
    - $`k_2 =  f(x^{n} + \frac{dt}{2} k_1)`$,
    - $`k_3 =  f(x^{n} + dt( 2k_2 - k_1))`$.

- Convergence order : 3.

## RK4 method

- Scheme:
$`x^{n+1} = x^{n} + \frac{dt}{6}  \left( k_1 + 2 k_2 + 2 k_3  + k_4\right)`$
  - with
    - $`k_1 =  f(x^{n+1})`$,
    - $`k_2 =  f(x^{n+1} + \frac{dt}{2} k_1)`$,
    - $`k_3 =  f(x^{n+1} + \frac{dt}{2} k_2)`$,
    - $`k_4 =  f(x^{n+1} + dt k_3)`$.

- Convergence order : 4.

