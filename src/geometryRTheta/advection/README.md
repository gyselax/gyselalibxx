# Advection operator

## Studied equation

The studied equation is the following 2D transport equation type :

```math
\partial_t f(t,x,y) + A(t,x,y)\cdot\nabla f(t,x,y) = 0,
```

with $`f(0,x,y) = f_0(x,y)`$ and *A* the advection field.

**We want to solve it on a polar grid so we have:**  $`(t,x,y) = (t,x(r,\theta),y(r,\theta))`$.

## Backward Semi-Lagrangian method

The method used to solve the equation is a backward Semi-Lagrangian method (BSL).
It uses the conservation property along the characteristics:

```math
\forall t, \quad f(t, x, y) = f(s, X(t; s, x, y), Y(t; s, x, y))
```

with:

```math
\begin{aligned}
\partial_t X (t; s, x, y) = A(t,X(t; s, x, y),Y(t; s, x, y)) \cdot e_x,\\
\partial_t Y (t; s, x, y) = A(t,X(t; s, x, y),Y(t; s, x, y)) \cdot e_y,\\
X(s; s, x, y) = x,\\
Y(s; s, x, y) = y.
\end{aligned}
```

So to compute the advected function at the next time step,

- we compute the feet of the characteristics $`X(t^n; t^{n+1}, x_i, y_j)`$ and $`Y(t^n; t^{n+1}, x_i, y_j)`$
 for each mesh point $`(x_i, y_j)`$ with a time integration method ;
- we interpolate the function $f(t = t^n)$ on the feet of the characteristics.
 The property ensures that the interpolation gives the function at the next time step $f(t = t^{n+1})$.

## Time integration methods

There are multiple time integration methods available which are implemented in the ITimeStepper child classes. For example:

- Explicit Euler method: Euler;
- Crank-Nicolson method: CrankNicolson;
- Runge-Kutta 3 method: RK3;
- Runge-Kutta 4 method: RK4;

We are listing the different schemes for this equation $`\partial_t X (t) = A(t, X(t))`$.

We write $X (t) = X (t; s, x)$, $X^n = X(t^n)$ and $A^n(X) = A(t^n, X)$ for a time discretisation
$`\{t^n; t^n > t^{n-1},  \forall n\}_n`$.

### Explicit Euler method

- Scheme:
$X^n = X^{n+1} - dt A^{n+1}(X^{n+1})$

- Convergence order : 1.

### Crank-Nicolson method

- Scheme:
$X^{k+1} = X^{n+1} - \frac{dt}{2} \left( A^{n+1}(X^{n+1}) + A^k(X^k) \right)$ and
$X^{n+1} = X^{k+1}$ once converged.

- Convergence order : 2.

### RK3 method

- Scheme:
$`X^n = X^{n+1} - \frac{dt}{6}  \left( k_1 + 4 k_2 + k_3 \right)`$
  - with
    - $`k_1 =  A^{n+1}(X^{n+1})`$,
    - $`k_2 =  A{n+1/2} (X^{n+1} - \frac{dt}{2} k_1)`$,
    - $`k_3 =  A{n+1/2} (X^{n+1} - dt( 2k_2 - k_1))`$.

- Convergence order : 3.

### RK4 method

- Scheme:
$`X^n = X^{n+1} - \frac{dt}{6}  \left( k_1 + 2 k_2 + 2 k_3  + k_4\right)`$
  - with
    - $`k_1 =  A^{n+1}(X^{n+1})`$,
    - $`k_2 =  A^{n+1/2} (X^{n+1} - \frac{dt}{2} k_1)`$,
    - $`k_3 =  A^{n+1/2} (X^{n+1} - \frac{dt}{2} k_2)`$,
    - $`k_4 =  A^{n} (X^{n+1} - dt k_3)`$.

- Convergence order : 4.

## Advection domain

There are two advection domains to consider:

- the physical domain;
- the pseudo-Cartesian domain.

It seems logical to use the **physical domain**, where the studied equation is given, as the advection domain.

However, we want to solve this equation on a polar grid. So before advecting, we have to
compute the mesh points in the physical domain using a mapping function $\mathcal{F}$:

```math
\mathcal{F} : (r,\theta)_{i,j} \mapsto  (x,y)_{i,j}.
```

This adds some steps to the advection operator, we now have to compute

- the mesh points in the physical domain using $\mathcal{F}$;
- the feet of the characteristics in the physical domain;
- the feet of the characteristics in the logical domain (polar grid) using $\mathcal{F}^{-1}$;
- then interpolate the advection function at the feet of the characteristics in the logical domain.

The third step can be difficult especially if the mapping function $\mathcal{F}$ is not analytically invertible.
It is not impossible, but the computations can be costly.

That is why, we introduce the **pseudo-Cartesian domain**.
We use another mapping function $\mathcal{G}$ such that:

```math
 \mathcal{G} : (r,\theta)_{i,j} \mapsto  (x,y)_{i,j} = (r\cos(\theta), r\sin(\theta))_{i,j}.
```

Then the four previous steps become

- calculate the mesh points in the pseudo-Cartesian domain using $\mathcal{G}$;
- calculate the advection field $A$ in the pseudo-Cartesian domain using the Jacobian matrix of $(\mathcal{F}\circ\mathcal{G}^{-1})^{-1}$;
- calculate the feet of the characteristics in the pseudo\_Cartesian domain;
- calculate the feet of the characteristics in the logical domain (polar grid) using $\mathcal{G}^{-1}$;
- interpolate the advection function at the feet of the characteristics in the logical domain.

Here, $\mathcal{G}$ is analytically invertible (we can fix  $\mathcal{G}^{-1}(x = 0, y = 0) = (r = 0, \theta = 0)$)
and  $`(J_{\mathcal{F}}J_{\mathcal{G}}^{-1})^{-1}`$ is well-defined. The details are given in Zoni et al. (2019) [^1].

**Remark 1:** if $\mathcal{F}$ is the circular mapping function, then the physical domain and the pseudo-Cartesian domain are the same.

**Remark 2:** if the mapping function is analytically invertible, it is less costly to advect in the physical domain.

## Advection Field

In the studied equation, the advection field is given along the physical domain axis:

```math
\partial_t f + A \cdot \nabla f = 0.
```

The BslAdvectionRTheta operator can take as input
the advection field expressed on the $`(e_x, e_y)`$ basis of the physical domain or
the advection field expressed on the $`(e_r, e_\theta)`$ [contravariant basis](../../../docs/standards/mathematical_and_physical_conventions.md) of the logical domain,

```math
A = A_x e_x + A_y e_y \quad \text{or} \quad A = A^r e_r + A^\theta e_\theta.
```

The advection field can be computed thanks to the AdvectionFieldFinder operator. This operator returns both expressions of the advection field (see [advection\_field\_rtheta](./../advection_field/README.md)).

- If the advection field is directly expressed on the $`(e_x, e_y)`$ basis of the physical domain, no treatment is needed in the BslAdvectionRTheta operator.

- If the advection field is expressed on the $`(e_r, e_\theta)`$ contravariant basis of the logical domain, then we need to compute the advection field on the $`(e_x, e_y)`$ basis to advect in the physical domain.

To pass from the composants on the $`(e_r, e_\theta)`$ contravariant basis to the composants on the $`(e_x, e_y)`$ basis, we use the Jacobian matrix *J* of the coordinate transformation $`(r,\theta) \mapsto (x,y)`$,

```math
\begin{bmatrix}
    A\cdot e_x \\
    A\cdot e_y
\end{bmatrix}
=
J
\begin{bmatrix}
    A \cdot e_{r} \\
    A \cdot e_{\theta}
\end{bmatrix}.
```

## Unit tests

The test of the advection operator are implemented in the `tests/geometryRTheta/advection_rtheta/` folder
([advection\_2d\_rp](./../../../tests/geometryRTheta/advection_rtheta/README.md)).

It tests:

- the 4 time integration methods;
- the different mappings and advection domains:
  - Circular mapping in the physical domain;
  - Czarny mapping in the physical domain;
  - Czarny mapping in the pseudo-Cartesian domain;
  - Discrete mapping of the Czarny mapping in the pseudo-Cartesian domain.
- on 3 different simulations:
  - simulation 1: translation of Gaussian function
    - $`f_0(x,y) = \exp\left( - \frac{(x- x_0)^2}{2 \sigma_x^2} - \frac{(y- y_0)^2}{2 \sigma_y^2} \right)`$,
    - $`A(t, x, y) = (v_x, v_y)`$ .
  - simulation 2: rotation of Gaussian function
    - $`f_0(x,y) = \exp\left( - \frac{(x- x_0)^2}{2 \sigma_x^2} - \frac{(y- y_0)^2}{2 \sigma_y^2} \right)`$,
    - $`A(t, x, y) = J_{\mathcal{F}_{\text{circular}}}(v_r, v_\theta)`$.
  - simulation 3: decentred rotation (test given in Zoni et al. (2019) [^1]).
    - $`f_0(x,y) = \frac{1}{2} \left( G(r_1(x,y)) + G(r_2(x,y))\right)`$,
      - with
        - $`G(r) = \cos\left(\frac{\pi r}{2 a}\right)^4 * 1_{r<a}(r)`$,
        - $`r_1(x, y) = \sqrt{(x-x_0)^2 + 8(y-y_0)^2}`$
        - $`r_2(x, y) = \sqrt{8(x-x_0)^2 + (y-y_0)^2}`$
    - $`A(t, x, y) = \omega(y_c - y, x - x_c)`$.

The tests of the convergence order are made for constant CFL which means it checks the slope of the errors
(infinity norm of the difference at the final time)
for $`(N_r\times N_\theta, dt) = (N_{r,0}\times N_{\theta,0}, dt_0)`$ and then $`(n*N_{r,0}\times n*N_{\theta,0}, dt_0/n)`$
for $n = 1, 2, 4, 8,  ...$.

## Contents

This folder contains:

- iadvection\_rtheta.hpp : define the base class for advection operator (IAdvectionRTheta).
  - bsl\_advection\_rtheta.hpp : define the advection operator described just before (BslAdvectionRTheta).

## References

[^1]: Edoardo Zoni, Yaman Güçlü. "Solving hyperbolic-elliptic problems on singular mapped
disk-like domains with the method of characteristics and spline finite elements".
([https://doi.org/10.1016/j.jcp.2019.108889](https://doi.org/10.1016/j.jcp.2019.108889).)
Journal of Computational Physics (2019).
