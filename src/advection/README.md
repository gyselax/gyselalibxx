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
DFieldX f(x_dom); 
DFieldX A(x_dom); 
using IDimInterest = IDimX; 
```


* Equation on a 2D domain (and time dimension): 
```math 
    \partial_t f(t,x,y) + A_x(x,y)\partial_{x}f(t,x,y) = 0,
```

Here $`\Omega' = \Omega \in \mathbb{R}^2`$. In the code, it would correspond to 
```cpp
DFieldXY f(xy_dom); 
DFieldXY A(xy_dom); 
using IDimInterest = IDimX; 
```


* Equation on a 2Dx2V domain (and time dimension): 
```math 
    \partial_t f(t,x,y,v_x,v_y) + A_x(x,y)\partial_{x}f(t,x,y,v_x,v_y) = 0,
```

Here $`\Omega' \in \mathbb{R}^2`$ and $`\Omega \in \mathbb{R}^4`$. In the code, it would correspond to 
```cpp
DFieldXYVxVy f(xyvxvy_dom); 
DFieldXY A(xy_dom); 
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
DFieldXVx f(xvx_dom); 
DFieldXVx A(xvx_dom); 
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
DFieldSpXVx f(sp_xvx_dom); 
DFieldSpX A(sp_x_dom); 
using IDimInterest = IDimX; 
```


### Parameters
The operator takes as templated parameters: 
* IDimInterest: a dimension of interest (or advection dimension) wich refers to the dimension along wich we advect the function; 
* AdvectionDomain: an advection domain, which refers to the domain where the advection field is defined. It has to contain the dimension of interest for the interpolation of the advection field in the time integration method; 
* FunctionDomain: the full domain where the function we want to advect is defined; 
* AdvectionFieldBuilder: a spline builder type for the advection field. 
* AdvectionFieldEvaluator: a spline evaluator type for the advection field. 
* TimeStepper: a time integration method (see [timestepper](./../timestepper/README.md)) to solve the characteristic equation. It has to be defined on the advection field domain. The feet have to be a Field of coordinates of the dimension of interest defined on the advection field domain.

**Remark/Warning:** the BslAdvection1D operator is built with builder and evaluator for the advection field and interpolator for the function we want to advect. Theses operators have to be defined on the same domain as the advection field and function. For instance, if the advection field and/or the function are defined on the species dimension, then the interpolators have to contain the species dimension in its batched dimensions (see tests in the `tests/advection/` folder). 

**Remark/Warning:** The advection field need to use interpolation on B-splines. So we cannot use other type of interpolator for the advection field. However there is no constraint on the interpolator of the advected function. 
