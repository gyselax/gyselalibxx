# Advection methods

The `advection/` folder gathers the backward semi lagrangian scheme classes. There are two main classes, AdvectionSpatial and AdvectionVelocity. Implementing these operators separately makes sense, since we are using time splitting.

The feet of the characteristic curves are used to interpolate the updated distribution function on mesh points. It uses batched interpolators so the interpolation step is done over the whole distribution function.

## Spatial advection
Here the purpose is the advection along a direction on the physical space dimension of the phase space.
The dynamics of the motion on the spatial dimension are governed by the following equation.

$$ \frac{df_s}{dt}= \sqrt{\frac{m_e}{m_s}} v \frac{\partial f_s}{\partial x} $$

## Velocity advection
Here the purpose is the advection along a direction on the velocity space dimension of the phase space.
The dynamics of the motion on the velocity dimension are governed by the following equation, where E is the electric field.

$$ \frac{df_s}{dt}= q_s \sqrt{\frac{m_e}{m_s}} E \frac{\partial f_s}{\partial v} $$


## 1D advection with a given advection field
The purpose of the BslAdvection1D operator is an advection along a given direction of the phase space. The advection field is given as input.
The dynamics of the motion are governed by the following equation.
```math
    \partial_t f(t,x) + A\partial_{x_i}(x')f(t,x) = 0,
    \qquad x \in \Omega, x' \in \Omega',
```

with
* $`f`$ the function to advect. It is defined on $`\Omega`$ domain and the time dimension;
* $`A`$ the advection field. It could be defined on a subdomain $`\Omega'\subset \Omega`$;
* and $`x_i`$ a given direction. The advection field domain has to be defined on this dimension for the time integration method.


### Example of use
Here are some examples of equation types the BslAdvection1D operator can solve:
* Equation on a 1D domain (and time dimension):
```math
    \partial_t f(t,x) + A_x(x)\partial_{x}f(t,x) = 0,
    \qquad x \in \Omega,
```

Here $`\Omega' = \Omega \in \mathbb{R}`$. In the code, it would correspond to
```cpp
DFieldX f(idx_range_x);
DFieldX A(idx_range_x);
using IDimInterest = IDimX;
```


* Equation on a 2D domain (and time dimension):
```math
    \partial_t f(t,x,y) + A_x(x,y)\partial_{x}f(t,x,y) = 0,
```

Here $`\Omega' = \Omega \in \mathbb{R}^2`$. In the code, it would correspond to
```cpp
DFieldXY f(idx_range_xy);
DFieldXY A(idx_range_xy);
using IDimInterest = IDimX;
```


* Equation on a 2Dx2V domain (and time dimension):
```math
    \partial_t f(t,x,y,v_x,v_y) + A_x(x,y)\partial_{x}f(t,x,y,v_x,v_y) = 0,
```

Here $`\Omega' \in \mathbb{R}^2`$ and $`\Omega \in \mathbb{R}^4`$. In the code, it would correspond to
```cpp
DFieldXYVxVy f(idx_range_xyvxvy);
DFieldXY A(idx_range_xy);
using IDimInterest = IDimX;
```


*  Equation on a 1Dx1V domain (and time dimension):
```math
\begin{aligned}
    & \partial_t f(t,x,v_x) + A_x(x,v_x)\partial_{x}f(t,x,v_x) = 0,
    \qquad \text{ with for instance, } A_x(x,v_x) = v_x, \\
    & \text{or } \partial_t f(t,x,v_x) + A_{v_x}(x,v_x)\partial_{v_x}f(t,x,v_x) = 0,
    \qquad \text{ with for instance, } A_{v_x}(x,v_x) = E(x),
\end{aligned}
```

Here $`\Omega' = \Omega \in \mathbb{R}^2`$. In the code, it would correspond to
```cpp
DFieldXVx f(idx_range_xvx);
DFieldXVx A(idx_range_xvx);
using IDimInterest = IDimX;
```
or
```cpp
using IDimInterest = IDimVx;
```

*  Equation on a 1Dx1V domain with species dimension (and time dimension):
```math
\begin{aligned}
    & \partial_t f_s(t,x,v_x) + A_{s,x}(x)\partial_{x}f_s(t,x,v_x) = 0,
\end{aligned}
```

Here $`\Omega' = \Omega \in \mathbb{R}^2`$. In the code, it would correspond to
```cpp
DFieldSpXVx f(idx_range_sp_xvx);
DFieldSpX A(idx_range_sp_x);
using IDimInterest = IDimX;
```


### Parameters
The operator takes as templated parameters:
* IDimInterest: a dimension of interest (or advection dimension) which refers to the dimension along which we advect the function;
* AdvectionDomain: an advection domain, which refers to the domain where the advection field is defined. It has to contain the dimension of interest for the interpolation of the advection field in the time integration method;
* FunctionDomain: the full domain where the function we want to advect is defined;
* AdvectionFieldBuilder: a spline builder type for the advection field.
* AdvectionFieldEvaluator: a spline evaluator type for the advection field.
* TimeStepper: a time integration method (see [timestepper](./../timestepper/README.md)) to solve the characteristic equation. It has to be defined on the advection field domain. The feet have to be a Field of coordinates of the dimension of interest defined on the advection field domain.

**Remark/Warning:** the BslAdvection1D operator is built with builder and evaluator for the advection field and interpolator for the function we want to advect. Theses operators have to be defined on the same domain as the advection field and function. For instance, if the advection field and/or the function are defined on the species dimension, then the interpolators have to contain the species dimension in its batched dimensions (see tests in the `tests/advection/` folder).

**Remark/Warning:** The advection field need to use interpolation on B-splines. So we cannot use other type of interpolator for the advection field. However there is no constraint on the interpolator of the advected function.

## PolarFootFinder

These methods are designed to calculate the foot of the characteristic on the polar plane for a 2D transport equation of the type:
```math
\partial_t f(t,x(r,\theta),y(r,\theta)) + A(t,x(r,\theta),y(r,\theta))\cdot\nabla f(t,x(r,\theta),y(r,\theta)) = 0,
```

with $`f(0,x,y) = f_0(x,y)`$ and *A* the advection field.

The characteristics are the solutions $X$ and $Y$ of the equations:
```math
\partial_t X (t; s, x, y) = A_x(t,X(t; s, x, y),Y(t; s, x, y)),\\
\partial_t Y (t; s, x, y) = A_y(t,X(t; s, x, y),Y(t; s, x, y)),\\
X(s; s, x, y) = x,\\
Y(s; s, x, y) = y.
```

The characteristic feet are calculated using a time integration method. For details of available methods see [timestepper](../timestepper/README.md).

### Advection domain

There are two advection domains to consider:
 - the physical domain;
 - the pseudo-Cartesian domain.

The logical domain is not used as it would be hard to calculate the feet close to the O-point in this domain.

It seems logical to use the **physical domain**, where the studied equation is given, as the advection domain.

However, we want to solve this equation on a polar grid. So before advecting, we have to
compute the mesh points in the physical domain using a mapping function $\mathcal{F}$:

```math
\mathcal{F} : (r,\theta)_{i,j} \mapsto  (x,y)_{i,j}.
```

This adds some steps to the advection operator, we now have to compute
 1. the mesh points in the physical domain using $\mathcal{F}$;
 2. the characteristic feet in the physical domain;
 3. the characteristic feet in the logical domain (polar grid) using $\mathcal{F}^{-1}$;
 4. then interpolate the advection function at the  characteristic feet in the logical domain.

The third step can be difficult especially if the mapping function $\mathcal{F}$ is not analytically invertible.
It is not impossible but the computations can be costly.


That is why, we introduce the **pseudo-Cartesian domain**.
We use another mapping function $\mathcal{G}$ such that:

```math
 \mathcal{G} : (r,\theta)_{i,j} \mapsto  (x,y)_{i,j} = (r\cos(\theta), r\sin(\theta))_{i,j}.
```

Then the four previous steps become
 1. calculate the mesh points in the pseudo-Cartesian domain using $\mathcal{G}$;
 2. calculate the advection field $A$ in the pseudo-Cartesian domain using the Jacobian matrix of $(\mathcal{F}\circ\mathcal{G}^{-1})^{-1}$;
 3. calculate the characteristic feet in the pseudo-Cartesian domain;
 4. calculate the characteristic feet in the logical domain (polar grid) using $\mathcal{G}^{-1}$;

Here, $\mathcal{G}$ is analytically invertible (we can fix  $`\mathcal{G}^{-1}(x = x_0, y = y_0) = (r = 0, \theta = 0)`$)
and  $`(J_{\mathcal{F}}J_{\mathcal{G}}^{-1})^{-1}`$ is well-defined. The details are given in Edoardo Zoni's article [1].


**Remark 1:** if $\mathcal{F}$ is the circular mapping function, then the physical domain and the pseudo-Cartesian domain are the same.

**Remark 2:** if the mapping function is analytically invertible, it is less costly to advect in the physical domain.

# References
[1] Edoardo Zoni, Yaman Güçlü. "Solving hyperbolic-elliptic problems on singular mapped
disk-like domains with the method of characteristics and spline finite elements".
([https://doi.org/10.1016/j.jcp.2019.108889](https://doi.org/10.1016/j.jcp.2019.108889).)
Journal of Computational Physics (2019).
