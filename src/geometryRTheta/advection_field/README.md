# Advection Field finder


The operator implemented here is a previous step to the advection operator. 
It computes the advection field along the axes in the physical domain or the axes in the logical domain.

Currently, the implemented case is:
* Guiding centre equations system.


## Guiding centre case

The studied equation system is of the following type : 
```math
\left\{
\begin{aligned}
\partial_t \rho + A\cdot\nabla \rho = 0, \\
A = E \wedge e_z, \\
E = - \nabla  \phi, \\
- \nabla \cdot \nabla \phi = \rho,
\end{aligned}
\right.
```

with $`\rho`$ the density, $`\phi`$ the electrostatic potential and $`E`$ the electrical field. 

The AdvectionFieldFinder computes the advection field $`A`$ from the electrical field $`\phi`$ returned by the PolarSplineFEMPoissonLikeSolver. 
It has two types of `operator()`: 
* one returning the advection field expressed on the vectors $`(e_x, e_y)`$ of the physical domain: $`A = A_x e_x + A_y e_y`$
* and another returning the advection field expressed on the vectors $`(e^r, e^\theta)`$ of the [contravariant basis](#docs_mathematical_and_physical_conventions) on the logical domain: $`A = A^r e_r + A^\theta e_\theta`$. 

The PolarSplineFEMPoissonLikeSolver can return the solution $`\phi`$ of the PDE under two forms:
* a Field of values of the solution on the mesh points of the grid; 
* a PolarSplineMem representation of the solution. 

The AdvectionFieldFinder can handle as input the two forms. 
If a Field is given as input, it computes the spline representation (on the cross-product of two 1D bases) using a SplineBuilder2D. 
The spline representation is needed to compute the derivatives of the function $`\phi`$. 
If the PolarSplineMem representation is given as input, it can directly compute the derivatives of the function $`\phi`$. 

Once the advection field computed, it is given as input to the BslAdvectionRTheta operator to advect the density $`\rho`$ function. 
The BslAdvectionRTheta operator can handle the advection with an advection field along $`(x,y)`$ and with an advection field along $`(r,\theta)`$. 
But as the BslAdvectionRTheta operator advects in the physical domain, it is recommended to work with the advection field along $`(x,y)`$.


### Advection field along the physical domain axis 

Thanks to the spline representation, the derivatives $`\partial_r \phi`$ and $`\partial_\theta \phi`$ are computed. 
The computation of the electrical field can be ill-defined around the O-point, so we treat this area separately. 

* If $`r > \varepsilon`$, we use 
```math
\nabla \phi
= \partial_x \phi e_x + \partial_y \phi e_y 
= J^{-T}(\partial_r \phi e^r + \partial_\theta \phi e^\theta), 
\qquad 
\text{in the covariant basis,}
```

with $`J`$  the Jacobian matrix and *G* the tensor metric of the mapping $`\mathcal{F}: (r,\theta)\mapsto(x,y)`$. 
Then the electric field is given by 
```math
E = -\nabla \phi
= - \partial_x \phi e_x - \partial_y \phi e_y 
```

and the advection field in the basis $`(e_x, e_y)`$ by
```math
A = E\wedge e_z 
= 
\begin{bmatrix}
     E\cdot e_y  \\
    -E\cdot e_x 
\end{bmatrix} 
= 
\begin{bmatrix}
    -\partial_y \phi \\
     \partial_x \phi
\end{bmatrix}. 
```

* If $`r \leq \varepsilon`$, we linearise. The method is detailed in [Zoni et al. (2019)](#zoni). We use only the derivatives along $`r`$ at two linearly independent directions of $`\theta`$ : $`\theta_1`$ and $`\theta_2`$
```math
\partial_r \phi (0, \theta_1) = \left[\partial_r x  \partial_x \phi + \partial_r y  \partial_y \phi \right] (0, \theta_1), \\
\partial_r \phi (0, \theta_2) = \left[\partial_r x  \partial_x \phi + \partial_r y  \partial_y \phi \right] (0, \theta_2).
```

From these equations, we deduce the (unique) values of $`\partial_x\phi`$ and $`\partial_y\phi`$ at $`(x,y) = (0,0)`$,

```math
\begin{bmatrix}
    \partial_x \phi (0, \theta) \\
    \partial_y \phi (0, \theta)
\end{bmatrix}
 = 
 \begin{bmatrix}
    \partial_r x (0, \theta_1)  & \partial_r y (0, \theta_1) \\
    \partial_r x (0, \theta_2)  & \partial_r y (0, \theta_2)
\end{bmatrix} ^{-1}
\begin{bmatrix}
    \partial_r \phi (0, \theta_1) \\
    \partial_r \phi (0, \theta_2)
\end{bmatrix}.
```

Then we compute $`E`$ at $`(x,y) = (0,0)`$ and $`(x,y) = \mathcal{F}(\varepsilon,\theta)`$ $`\forall \theta`$ (for $`\varepsilon\neq 0`$, we use the Jacobian matrix as previously) and we linearise

```math
E(r, \theta) = \left( 1 - \frac{r}{\varepsilon} \right)  E(0, \theta) + \frac{r}{\varepsilon} E(\varepsilon, \theta), \\
```

As previously, we compute the advection field by 
```math
A = E\wedge e_z 
= 
\begin{bmatrix}
      E\cdot e_y  \\
    - E\cdot e_x 
\end{bmatrix}.
```

(In the code, we chose $`\theta_1 = \frac{\pi}{4}`$ and $`\theta_2  = - \frac{\pi}{4}`$, and $\varepsilon = 10^{-12}$.)


### Advection field along the logical domain axis

Firstly, the derivatives $`\partial_r \phi`$ and $`\partial_\theta \phi`$ are also computed here. 

#### General coordinates system 
* In a **general coordinate system**, the gradient of a scalar function in the logical domain is given in the [covariant basis](#docs_mathematical_and_physical_conventions) by 

```math
\nabla f = \sum_i \partial_{q_i} f b^j
```

with 
* $`J`$ the Jacobian matrix associated with the mapping function of the system $`\mathcal{F}:(q_1,..., q_N)\mapsto(x_1, ..., x_N)`$, 
* $`G = J^T J`$ the metric tensor, whose components are $`g_{ij}`$,  
* $`G^{-1}`$ the inverse metric tensor, whose components are $`g^{ij}`$, 
* and $`b^j`$ the unnormalised local covariant vectors. 

In 2D with $`(q_1, q_2) = (r,\theta)`$ and $`(x_1, x_2) = (x,y)`$, it can be rewritten as the following matrix system 
```math
\hat{\nabla} f 
= \partial_x f e_x + \partial_y f e_y
= \partial_r f e^r + \partial_\theta f e^\theta
```

With the composants linked by the following relation, 
```math
\begin{bmatrix}
    \partial_{x} f \\
    \partial_{y} f
\end{bmatrix}
= 
J^{-T}
\begin{bmatrix}
    \partial_{r} f \\
    \partial_{\theta} f
\end{bmatrix}.
```


#### Application to the advection field
The gradient is expressed in the *covariant* basis. We express the advection field in the *contravariant* basis to use the nice property $`e_r\cdot e^r = 1`$ and $`e_\theta\cdot e^\theta = 1`$. 
* In the [contravariant basis](#docs_mathematical_and_physical_coneventions) $`(e_r, e_\theta)`$, 
we compute the electric field,
```math
E
= 
\begin{bmatrix}
    E \cdot e_r \\
    E \cdot e_{\theta}
\end{bmatrix}
= 
G^{-1}
\begin{bmatrix}
    - \partial_{r} \phi \\
    - \partial_{\theta} \phi
\end{bmatrix}.
```

with $`G^{-1}`$ the inverse [metric tensor](#docs_mathematical_and_physical_coneventions__Metric_tensor). 

Then the advection field is given by 
```math
A = E \wedge e_z = 
\begin{bmatrix}
     E \cdot e_y \\
    -E \cdot e_x
\end{bmatrix},
```

with in contravariant basis, 
```math
\begin{bmatrix}
    E \cdot e_x \\
    E \cdot e_y
\end{bmatrix}
= 
J
\begin{bmatrix}
    E \cdot e_r \\
    E \cdot e_\theta
\end{bmatrix} ,
\qquad
\begin{bmatrix}
    A \cdot e_x \\
    A \cdot e_y
\end{bmatrix}
= 
J
\begin{bmatrix}
    A \cdot e_r \\
    A \cdot e_\theta
\end{bmatrix}. 
```

So, 
```math
\begin{bmatrix}
    A \cdot e_r \\
    A \cdot e_\theta
\end{bmatrix}
= J^{-1}
\begin{bmatrix}
    A \cdot e_x \\
    A \cdot e_y
\end{bmatrix}
= J^{-1}
\begin{bmatrix}
     E \cdot e_y \\
    -E \cdot e_x
\end{bmatrix}
= J^{-1}
\begin{bmatrix}
    J_{2,1} & J_{2,2} \\
    -J_{1,1} & -J_{1,2}
\end{bmatrix}
\begin{bmatrix}
    E \cdot e_r \\
    E \cdot e_\theta
\end{bmatrix}
= \frac{1}{\det(J)}
\begin{bmatrix}
    (J_{11}J_{12} + J_{21}J_{22})  &  (J_{22}^2 + J_{12}^2) \\
    -(J_{11}^2 + J_{21}^2)         & -(J_{11}J_{12} + J_{21}J_{22})
\end{bmatrix}
\begin{bmatrix}
    E \cdot e_r \\
    E \cdot e_\theta
\end{bmatrix}
```

with $`J_{ij}`$ the elements of the matrix *J*. 

Warning, the matrix $`G^{-1}`$ can be ill-defined for $r = 0$. 

*Example: circular mapping:* 
```math
G^{-1}
= 
\begin{bmatrix}
    1 & 0 \\
    0 & \frac{1}{r^2}
\end{bmatrix}.
```

In the code, the O-point is differently treated. The domain is split between a domain without the O-point ($`(0,\theta), \forall \theta`$) and the domain containing only the O-point. For the first domain, we compute the advection field along the logical axis as explain previously. On the second domain, we compute the unique value of the advection field along the physical axis using the linearisation done in the [Advection field along the physical domain axis](#src_geometryRTheta_advection_field__Guiding_centre_case) section.



# References 

<a name="zoni"></a> [1] Edoardo Zoni, Yaman Güçlü, "Solving hyperbolic-elliptic problems on singular mapped disk-like domains with the 
method of characteristics and spline finite elements", https://doi.org/10.1016/j.jcp.2019.108889, Journal of Computational Physics, 2019.


# Contents

* advection\_field\_rtheta.hpp : containing AdvectionFieldFinder with the advection field computation for the guiding centre simulation. 
