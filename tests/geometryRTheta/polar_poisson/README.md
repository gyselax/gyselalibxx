# Tests on the 2D polar poisson solver


The tests implemented in this folder test the 2D polar poisson solver implemented in the `src/geometryRTheta/poisson/` folder 
( [poisson](./../../../src/geometryRTheta/poisson/README.md) ).


## Polar Poisson solver

The PolarSplineFEMPoissonLikeSolver is tested on a circular mapping (CircularToCartesian) and on a Czarny mapping (CzarnyToCartesian) with
(the test cases are given in Emily Bourne's thesis [1]). 

The studied equation is

$$L\phi = - \nabla \cdot (\alpha \nabla \phi) + \beta \phi = \rho, $$

with  
 * $`\alpha(r) = \exp\left( - \tanh\left( \frac{r - r_p}{\delta_r} \right) \right)`$, with $`r_p = 0.7`$ and $`\delta_r = 0.05`$, 
 * $\beta(r) = \frac{1}{\alpha(r)}$. 


The defined test functions are: 

 * Polar solution: $\phi(x, y) = C r(x,y)^6 (r(x,y) -1)^6 \cos(m\theta)$,
 	* with $C = 2^{12}1e-4$ and $m = 11$; 
 * Cartesian solution: $\phi(x,y) = C (1+r(x,y))^6  (1 - r(x,y))^6 \cos(2\pi x) \sin(2\pi y)$, 
 	* with  $C = 2^{12}1e-4$. 
 	
 	
The VlasovPoissonSolver is also tested on a circular mapping (CircularToCartesian) and on a Czarny mapping (CzarnyToCartesian) 
for the same Poisson equation with the Cartesian solution. 
 		
 		
 		
 		
 		

## References 

[1] Emily Bourne, "Non-Uniform Numerical Schemes for the Modelling of Turbulence in the 5D GYSELA Code". December 2022.



## Contents

 * polarpoissonfemsolver.cpp : it tests the PolarSplineFEMPoissonLikeSolver. It solves the Poisson equation for a selected test case.
 * test\_cases.hpp : it defines RHS (ManufacturedPoissonTest) and exact solutions (PoissonSolution) of the Poisson equation. 
 * vlasovpoissonsolver.cpp : it tests the VlasovPoissonSolver. It solves the Poisson equation for a selected test case and computes the electric field in the physical domain.
 * poisson.yaml : the parameters of the tests.
 * test\_poisson.py : it launches twice polarpoissonfemsolver.cpp test and checks the convergence order of the solution.  
 * test\_vlasov\_poisson.py : it launches twice vlasovpoissonsolver.cpp test and checks the convergence order or the solution and the electric field. 


