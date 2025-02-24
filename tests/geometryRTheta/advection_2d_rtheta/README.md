# Tests on the 2D polar advection operator
 
## Tests on the 2D polar advection operator

The tests implemented in this folder test the 2D polar advection operator implemented in the `src/geometryRTheta/advection/` folder 
( [advection](./../../../src/geometryRTheta/advection/README.md) ).

The tests are made for different parameters which are:

### - The mapping and the domain used for the advection: 
 - Circular mapping in the physical domain (CircularToCartesian and AdvectionPhysicalDomain); 
 - Czarny mapping in the physical domain (CzarnyToCartesian and AdvectionPhysicalDomain); 
 - Czarny mapping in the pseudo-Cartesian domain (CzarnyToCartesian and AdvectionPseudoCartesianDomain); 
 - Discrete mapping of the Czarny mapping in the pseudo-Cartesian domain (DiscreteToCartesian and AdvectionPseudoCartesianDomain).
 	
### - The time integration method used to solve the characteristic equation: 
 - Explicit Euler (Euler); 
 - Crank-Nicolson (CrankNicolson); 
 - Runge-Kutta 3 (RK3); 
 - Runge-Kutta 4 (RK4). 
 	
### - The test simulation: 
 - simulation 1: translation of Gaussian function (TranslationSimulation)
   - $`f_0(x,y) = \exp\left( - \frac{(x- x_0)^2}{2 \sigma_x^2} - \frac{(y- y_0)^2}{2 \sigma_y^2} \right)`$, 
   - $`V(t, x, y) = (v_x, v_y)`$ . 
 - simulation 2: rotation of Gaussian function (RotationSimulation)
  - $`f_0(x,y) = \exp\left( - \frac{(x- x_0)^2}{2 \sigma_x^2} - \frac{(y- y_0)^2}{2 \sigma_y^2} \right)`$, 
  - $`V(t, x, y) = J_{\mathcal{F}_{\text{circular}}}(v_r, v_\theta)`$. 
 - simulation 3: decentred rotation (test given in Edoardo Zoni's article [1]) (DecentredRotationSimulation)
  - $`f_0(x,y) = \frac{1}{2} \left( G(r_1(x,y)) + G(r_2(x,y))\right)`$,
  - with 
     - $`G(r) = \cos\left(\frac{\pi r}{2 a}\right)^4 * 1_{r<a}(r)`$, 
     - $`r_1(x, y) = \sqrt{(x-x_0)^2 + 8(y-y_0)^2}`$ 
     - $`r_2(x, y) = \sqrt{8(x-x_0)^2 + (y-y_0)^2}`$ 
  - $`V(t, x, y) = \omega(y_c - y, x - x_c)`$. 


## Python tests

- animated\_curves.py: create `.mp4` video(s) of the advected function for the selected configuations among the 48 test ones.
	- Command to launch the test in this folder: `python3 animated_curves.py ../../../build/tests/geometryRTheta/advection_2d_rtheta/advection_ALL`
	or `python3 animated_curves.py ../../../build/tests/geometryRTheta/advection_2d_rtheta/<selected advection test case>`

- display\_all\_errors\_for\_gtest.py: Google test which tests the convergence order for the 48 configurations.
	- Command to launch the test in this folder: `python3 display_all_errors_for_gtest.py ../../../build/tests/geometryRTheta/advection_2d_rtheta/advection_ALL`

- display\_curves.py: display the curve of the function for 9 time steps between the initial and the final state.  
	- Command to launch the test in this folder: `python3 display_curves.py ../../../build/tests/geometryRTheta/advection_2d_rtheta/<selected advection test case>`

- display\_feet\_errors.py: compute the convergence order of the characteristic feet and display the computed and the exact feet. 
	- Command to launch the test in this folder: `python3 display_feet_errors.py ../../../build/tests/geometryRTheta/advection_2d_rtheta/<selected advection test case>`

- advection\_functions.py: define all the useful functions used in the other python files. 



## References
[1] Edoardo Zoni, Yaman Güçlü. "Solving hyperbolic-elliptic problems on singular mapped 
disk-like domains with the method of characteristics and spline finite elements". 
([https://doi.org/10.1016/j.jcp.2019.108889](https://doi.org/10.1016/j.jcp.2019.108889).)
Journal of Computational Physics (2019).

## Contents

- advection\_all\_tests.cpp : launch the 48 configurations on an uniform mesh. 
- advection\_selected\_test.cpp : launch 1 configuration. The selection is made in the CMakeLists.txt file.
- advection\_maths\_tools.cpp : define functions used in the test files. 
- test\_cases.hpp : define the three test simulations. 


