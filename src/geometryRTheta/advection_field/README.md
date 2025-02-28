# Advection Field finder


The operator implemented here is a previous step to the advection operator. 
It computes the advection field along the axes in the physical domain or the axes in the logical domain.

Currently, the implemented case is:
* Guiding centre equations system.


## Guiding centre case

The studied equation system is of the following type : 
```math
\partial_t \rho + A\cdot\nabla \rho = 0, \\
A = E \wedge e_z, \\
E = - \nabla  \phi, \\
- \nabla \cdot \nabla \phi = \rho,
```

with $`\rho`$ the density, $`\phi`$ the electrostatic potential and $`E`$ the electrical field. 

The AdvectionFieldFinder computes the advection field $`A`$ from the electrical field $`\phi`$ returned by the PolarSplineFEMPoissonLikeSolver. 
It has two types of `operator()`: 
* one returning the advection field along the axis of the physical domain: $`A = (A_x, A_y)`$
* and the another returning the advection field along the axis of the logical domain: $`A = (A_r, A_\theta)`$. 

The PolarSplineFEMPoissonLikeSolver can return the solution $`\phi`$ of the PDE under two forms:
* a Field of values of the solution on the mesh points of the grid; 
* a PolarSplineMem representation of the solution. 

The AdvectionFieldFinder can handle as input the two forms. 
If a Field is given as input, it computes the spline representation (on the cross-product of two 1D bases) using a SplineBuilder2D. 
The spline representation is needed to compute the derivatives of the function $`\phi`$. 
If the PolarSplineMem representation is given as input, it can directly compute the derivatives of the function $`\phi`$. 

Once the advection field computed, it is given as input to the BslAdvectionRP operator to advect the density $`\rho`$ function. 
The BslAdvectionRP operator can handle the advection with an advection field along $`(x,y)`$ and with an advection field along $`(r,\theta)`$. 
But as the BslAdvectionRP operator advects in the physical domain, it is recommend to work with the advection field along $`(x,y)`$.


### Advection field along the physical domain axis 

Thanks to the spline representation, the derivatives $`\partial_r \phi`$ and $`\partial_\theta \phi`$ are computed. 
The computation of the electrical field can be ill-defined around the O-point so we treat this area separately. 

* If $`r > \varepsilon`$, we use 
```math
\begin{bmatrix}
    \partial_x \phi \\
    \partial_y \phi \\
\end{bmatrix} 
= 
J^{-T}
\begin{bmatrix}
    \partial_r \phi \\
    \partial_\theta \phi \\
\end{bmatrix}
```

with $`J`$  the Jacobian matrix of the mapping $`\mathcal{F}: (r,\theta)\mapsto(x,y)`$. Then the electric field is given by 
```math
E = -\nabla \phi
= 
\begin{bmatrix}
    - \partial_x \phi \\
    - \partial_y \phi \\
\end{bmatrix} 
```

and the advection field by 
```math
A = E\wedge e_z 
= 
\begin{bmatrix}
    - E_y  \\
    E_x  \\
\end{bmatrix} 
= 
\begin{bmatrix}
    \partial_y \phi \\
    - \partial_x \phi \\
\end{bmatrix}. 
```

* If $`r \leq \varepsilon`$, we linearise. The method is detailed in Edoardo Zoni's article [1]. We use only the derivatives along $`r`$ at two  linearly independent directions of $`\theta`$ : $`\theta_1`$ and $`\theta_2`$
```math
\partial_r \phi (0, \theta_1) = \left[\partial_r x  \partial_x \phi + \partial_r y  \partial_y \phi \right] (0, \theta_1), \\
\partial_r \phi (0, \theta_2) = \left[\partial_r x  \partial_x \phi + \partial_r y  \partial_y \phi \right] (0, \theta_2).
```

From these equations, we deduce the (unique) values of $`\partial_x\phi`$ and $`\partial_y\phi`$ at $`(x,y) = (0,0)`$,

```math
\begin{bmatrix}
    \partial_x \phi (0, \theta) \\
    \partial_y \phi (0, \theta) \\
\end{bmatrix}
 = 
 \begin{bmatrix}
    \partial_r x (0, \theta_1)  & \partial_r y (0, \theta_1) \\
    \partial_r x (0, \theta_2)  & \partial_r y (0, \theta_2) \\
\end{bmatrix} ^{-1}
\begin{bmatrix}
    \partial_r \phi (0, \theta_1)  \\
   \partial_r \phi (0, \theta_2) \\
\end{bmatrix}.
```

Then we compute $`E`$ at $`(x,y) = (0,0)`$ and $`(x,y) = \mathcal{F}(\varepsilon,\theta)`$ $`\forall \theta`$ (for $`\varepsilon\neq 0`$, we use the Jacobian matrix as previously) and we linearise

```math
E_x(r, \theta) = \left( 1 - \frac{r}{\varepsilon} \right)  E_x(0, \theta) + \frac{r}{\varepsilon} E_x(\varepsilon, \theta), \\
E_y(r, \theta) = \left( 1 - \frac{r}{\varepsilon} \right)  E_y(0, \theta) + \frac{r}{\varepsilon} E_y(\varepsilon, \theta), 
```

As previously, we compute the advection field by 
```math
A = E\wedge e_z 
= 
\begin{bmatrix}
    - E_y  \\
    E_x  \\
\end{bmatrix} 
= 
\begin{bmatrix}
    \partial_y \phi \\
    - \partial_x \phi \\
\end{bmatrix}. 
```

(In the code, we chose $`\theta_1 = \frac{\pi}{4}`$ and $`\theta_2  = - \frac{\pi}{4}`$, and $\varepsilon = 10^{-12}$.)


### Advection field along the logical domain axis

Firstly, the derivatives $`\partial_r \phi`$ and $`\partial_\theta \phi`$ are also computed here. 

#### General coordinates system 
* In **general coordinates system**, the gradient of a function is given by 

```math
\nabla f = \sum_i \sum_j \partial_{x_i} f g^{ij} \sqrt{g_{jj}} \hat{e}_j, 
```

with 
* $`J`$ the Jacobian matrix associated the the mapping function of the system $`\mathcal{F}:(x_1, x_2)\mapsto(y_1,y_2)`$, 
* $`G = J^T J = [g_{ij}]_{ij}`$ the metric tensor, 
* $`G^{-1} = [g^{ij}]_{ij}`$ the inverse metric tensor 
* and $`\hat{e}_j`$ the normalised covariant vectors. 

In 2D, it can be rewritten as the following matrix system 
```math
\nabla f = 
D_{G} G^{-T}
\begin{bmatrix}
    \partial_{x_1} f \\
    \partial_{x_2} f \\
\end{bmatrix}
= 
\begin{bmatrix}
    \sqrt{g_{11}} & 0 \\
    0 & \sqrt{g_{22}} \\
\end{bmatrix}
\begin{bmatrix}
    g^{11} & g^{21} \\
    g^{12} & g^{22} \\
\end{bmatrix}
\begin{bmatrix}
    \partial_{x_1} f \\
    \partial_{x_2} f \\
\end{bmatrix}.
```

**Remark:** We can prove that $`\det(D_{G} G^{-T}) = 0`$ if $`\det(J) = 0`$ (if the coefficients of the matrix $`D_G`$ are not null).
So for an invertible matrix, we also have the relation 
```math
\begin{bmatrix}
    \partial_{x_1} f \\
    \partial_{x_2} f \\
\end{bmatrix}
= 
G^{T}D_{G}^{-1} \nabla f. 
```

From the relation 
```math
\begin{bmatrix}
    \partial_{y_1} f \\
    \partial_{y_2} f \\
\end{bmatrix}
= 
J^{-T}
\begin{bmatrix}
    \partial_{x_1} f \\
    \partial_{x_2} f \\
\end{bmatrix}, 
```

we deduce the following relation for invertible case
```math
\nabla_{y_1, y_2} f
= 
J D_{G}^{-1}
\nabla_{x_1, x_2} f,
```

with $`\nabla_{y_1, y_2} f = [\partial_{y_1} f, \partial_{y_2} f]^T`$ and $`\nabla_{x_1, x_2} f = \sum_i \sum_j \partial_{x_i} f g^{ij} \sqrt{g_{jj}} \hat{e}_j`$.


#### Application to the advection field
* In our case, we use this formula to compute the electric field along the logical axis: 
```math
E
= -\nabla \phi  
= - D_{G} (G^{-1})^{T}
\begin{bmatrix}
    \partial_{r} \phi \\
    \partial_{\theta} \phi \\
\end{bmatrix}.
```

Then the advection field is given by 
```math
A
= E \wedge e_z
= 
\begin{bmatrix}
    -E_{\theta} \\
    E_r \\
\end{bmatrix}.
```

Warning, the matrix $`(G^{-1})^{T}`$ is ill-defined for $r = 0$. 

*Example: circular mapping:* 
```math
(G^{-1})^{T}
= 
\begin{bmatrix}
    1 & 0 \\
    0 & \frac{1}{r^2} \\
\end{bmatrix}.
```

In the code, the O-point is differently treated. The domain is split between a domain without the O-point ($`(0,\theta), \forall \theta`$) and the domain containing only the O-point. For the first domain, we compute the advection field along the logical axis as explain previously. On the second domain, we compute the unique value of the advection field along the physical axis using the linearisation done in the *Advection field along the physical domain axis* section. 



# References 

[1] Edoardo Zoni, Yaman Güçlü, "Solving hyperbolic-elliptic problems on singular mapped disk-like domains with the 
method of characteristics and spline finite elements", https://doi.org/10.1016/j.jcp.2019.108889, Journal of Computational Physics, 2019.


# Contents

* advection\_field\_rp.hpp : containing AdvectionFieldFinder with the advection field computation for the guiding centre simulation. 
