# Polar Poisson solver


## The Poisson equation 

(For more details, see Emily Bourne's thesis "Non-Uniform Numerical Schemes for the Modelling of Turbulence
in the 5D GYSELA Code". December 2022.)

The equation we are solving here is 

```math
L\phi = - \nabla \cdot (\alpha \nabla \phi) + \beta \phi = \rho
```

with the boundary condition $\phi = 0$, on  $\partial \Omega$. 


To solve this equation, the PolarSplineFEMPoissonSolver uses a finite element method on the B-splines. 

### B-splines FEM

#### The B-splines 
We introduce a basis of B-splines $`\{B_l\}_l`$, cross-product of two 1D B-splines basis: 

```math
\{B_l\}_l = \{b_{i,r} b_{j,\theta}\}_{i n_{n,\theta} +j}.
```


#### Treatment of the O-point
 
To conserve $\mathcal{C}^k$ property at the O-point, the B-splines basis is modified. 

The $\frac{(k+1)(k+2)}{2}$ first elements in the 1D B-splines basis in the $r$ dimension are removed, and 
replaced by 3 functions, linear combinations of these removed B-splines. 
(See (2.34) in Emily Bourne's thesis.)

The linear composition uses as coefficients, the Berstein polynomials defined from barycentric coordinates. 

The B-splines with the treatment of the O-point, are called PolarBSplines and are written $`\{\hat{B}_l\}_l`$. 


#### Weak formulation
The Poisson equation is solved by solving its weak form: 

```math
\int_{\Omega} \lbrack \beta(r) \phi(r,\theta) \hat{B}_l(r,\theta) + \alpha(r) \nabla \phi(r,\theta) \cdot  \nabla \hat{B}_l(r,\theta) \rbrack |det(J_{\mathcal{F}}(r,\theta))| dr d\theta =  \int_{\Omega} \rho(r,\theta) \hat{B}_l(r,\theta) |det(J_{\mathcal{F}}(r,\theta))| dr d\theta
```

with $`\{\hat{B}_l\}_l`$ the B-splines basis. 


Written as matrix form, it gives

```math
(M + S) \hat{\phi} = M \hat{\rho},
```

with 
 * the mass matrix, $`M_{i,j} =  \int_{\Omega} \beta(r,\theta) \hat{B}_i(r,\theta)\hat{B}_j(r,\theta) |det(J_{\mathcal{F}}(r,\theta))| dr d\theta`$, 
 * the stiffness matrix, $`S_{i,j} =  \int_{\Omega} \alpha(r) \left[\sum_{\xi_1\in[r,\theta]} \sum_{\xi_2\in[r,\theta]} \partial_{\xi_1}\hat{B}_i(r,\theta) \partial_{\xi_2}\hat{B}_j(r,\theta) g^{\xi_1, \xi_2}\right] |det(J_{\mathcal{F}}(r,\theta))| dr d\theta`$
 	* with $`g^{\xi_1, \xi_2}`$, the scalar product between the unit vector in $`\xi_1`$ and the unit vector in $`\xi_2`$. 
 * the solution vector, $`\hat{\phi}_i = \sum_{l = 0}^{n_{b,\theta}(n_{b,r} -2) +3} \phi_l \hat{B}_l(r_i,\theta_i)`$, 
 * and the rhs vector,  $`\hat{\rho}_j = \sum_{l = 0}^{n_{b,\theta}(n_{b,r} -2) +3} \rho_l \hat{B}_l(r_j,\theta_j)`$. 

So we compute the solution B-splines coefficients $`\{\phi_l\}_l`$ by solving this matrix equation.  





## Unit tests 

The test are implemented in the `tests/geometryRTheta/polar_poisson/` folder 
([polar\_poisson](./../../../tests/geometryRTheta/polar_poisson/README.md)).

The PolarSplineFEMPoissonSolver is tested on a circular mapping (CircularToCartesian) and on a Czarny mapping (CzarnyToCartesian) with
(the test cases are given in Emily Bourne's thesis [1])
 
 * the Poisson coefficients: 
 	* $`\alpha(r) = \exp\left( - \tanh\left( \frac{r - r_p}{\delta_r} \right) \right)`$, with $`r_p = 0.7`$ and $`\delta_r = 0.05`$, 
 	* $\beta(r) = \frac{1}{\alpha(r)}$. 
 * for the following test functions: 
 	* polar solution: $\phi(x, y) = C r(x,y)^6 (r(x,y) -1)^6 \cos(m\theta)$,
 		* with $C = 2^{12}1e-4$ and $m = 11$; 
 	* cartesian solution: $\phi(x,y) = C (1+r(x,y))^6  (1 - r(x,y))^6 \cos(2\pi x) \sin(2\pi y)$, 
 		* with  $C = 2^{12}1e-4$. 
 		


## References 

[1] Emily Bourne, "Non-Uniform Numerical Schemes for the Modelling of Turbulence in the 5D GYSELA Code". December 2022.

[2] Edoardo Zoni, Yaman Güçlü, "Solving hyperbolic-elliptic problems on singular mapped disk-like domains with the 
method of characteristics and spline finite elements", https://doi.org/10.1016/j.jcp.2019.108889, Journal of Computational Physics, 2019.


## Contents

 * ipoissonsolver.hpp : Define a base class for the Poisson solvers: IPoissonSolver.
 * polarpoissonsolver.hpp : Define a Poisson solver using FEM on B-splines: PolarSplineFEMPoissonSolver. 
 * poisson\_rhs\_function.hpp : Define a rhs object (PoissonRHSFunction) for the Poisson equation (mainly used for vlasovpoissonsolver.hpp): PoissonRHSFunction. 
