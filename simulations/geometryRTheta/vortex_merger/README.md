# Vortex merger

## Studied problem
### Field setup 
We work with cylindrical coordinates $(r,\theta,z)$ on the corresponding orthonormal basis vectors $`(e_r, e_\theta, e_z)`$.


### General case
Conservation of the charge density equation and Maxwell-Poisson equation:

```math
\partial_t \rho + div(\rho v) = 0
```
```math
\Delta \phi + \frac{\rho}{\varepsilon_0} = 0
```
```math
v = \frac{E\wedge B}{|B|^2}
```
```math
E = -\nabla \phi
```

where
* $\rho = \rho(t,r,\theta)$ is the **electric charge density** of the particles (electrons with mass $`m_e`$ and charge $q$);
* $v$ is the velocity of the particles ($div(v) = 0$); 
* $\phi$ is the usual electric potential associated with the total electric field; 
* $(E,B)$ the electromagnetic field such that $`B = B_z(r) e_z`$ the external magnetic field


### Simplified case
We study a simplified case where $`B = e_z`$, $`\varepsilon_0 = 1`$. 
On a circular mapping, it gives 
$`v = -\nabla\phi \wedge e_z`$, $div(v) = 0$  
and with $`(e_r, e_\theta)`$ 
the unnormalized local contravariant base, 
```math
\partial_t \rho - \frac{\partial_\theta\phi}{r}\partial_r\rho + \partial_r\phi \frac{\partial_\theta\rho}{r} =0,
```

```math
-\Delta \phi = \rho.
```


## Test case - vortex merger
The test case implemented tests a vortex merger (see [1] article).

We suppose as perturbed initial condition

```math
\rho(0,r, \theta) = \rho_{eq}(x,y) + \varepsilon 
	\left( 
		\exp\left[ - \frac{(x - x_1^*)^2 + (y - y_1^*)^2}{2 \sigma^2} \right]
		+ \exp\left[ - \frac{(x - x_2^*)^2 + (y - y_2^*)^2}{2 \sigma^2} \right]
	\right)
```


with 
* $`(x_1^*, y_1^*) = (0.08, -0.14)`$, 
* $`(x_2^*, y_2^*) = (-0.08, 0.14)`$,
* $\sigma = 0.08$ and 
* $\varepsilon = 10^{-4}$.

The equilibrium solution is computed during the initialisation phase in the VortexMergerEquilibria class. 
The details are documented in [initialisation](./../../../src/geometryRTheta/initialisation/README.md). 


## References
[1]    Edoardo Zoni, Yaman Güçlü, "Solving hyperbolic-elliptic problems on singular mapped disk-like domains with the 
method of characteristics and spline finite elements", https://doi.org/10.1016/j.jcp.2019.108889, Journal of Computational Physics, 2019.

## Useful references 
[2]    Emily Bourne. “Non-Uniform Numerical Schemes for the Modelling of Turbulence in the 5D GYSELA Code”. PhD thesis. Aix-Marseille Université, Dec. 2022.

[3]    Eric Sonnendrucker, Jean Roche, Pierre Bertrand and Alain Gizzo. “The Semi-Lagrangian Method for the Numerical Resolution of the Vlasov Equation”. Journal of Computational Physics (1999).






## Contents

 - vortex\_merger.cpp : runs a vortex-merger simulation.
 - params.yaml : defines the parameters of the simulation. 
 
 Python files: (in `post_process/PythonScripts/geometryRTheta` folder)
 - animation\_rho\_phi.py : plot and save the density and electrical potential function in time. 
 - mass\_conservation.py : plot the relative errors of the mass of the particles. 
 
 
 ### Recommended parameters: 
* `Mesh:`
  * `r_size: 128` : number of cells in $r$-dimension. (Tests in Edoardo Zoni's article.)
  * `theta_size: 256` : number of cells in $\theta$-dimension. (Tests in Edoardo Zoni's article.)
  * `r_min: 0.0`  : start of $`r`$ domain. (Tests in Edoardo Zoni's article.)
  * `r_max: 1.0` : end of $`r`$ domain. (Tests in Edoardo Zoni's article.)

* `Time:`
  * `delta_t: 0.1` : time step. (Tests in Edoardo Zoni's article.)
  * `final_T: 10.0`: final time of the simulation. (Tests in Edoardo Zoni's article.)
  
* `Perturbation:`
  * `eps: 0.0001` : amplitude of the perturbation. 
  * `x_star_1: 0.08` : $`x`$ coordinate of the centre of the initial first vortex.
  * `y_star_1: -0.14` : $`y`$ coordinate of the centre of the initial first vortex.
  * `x_star_2: -0.08` : $`x`$ coordinate of the centre of the initial second vortex.
  * `y_star_2: 0.14` : $`y`$ coordinate of the centre of the initial second vortex.
  * `sigma: 0.08` : the standard deviation of the initial Gaussian function. 
  
* `Output:`
  * `time_step_diag: 5` : number of time steps between two recordings of the data. 
 
 
 ### Executables:
* `vortex_merger` : the vortex-merger simulation with the given parameters. 
 
The results of the simulation are saved in a `output` folder that the executable creates. 

The path to the `params.yaml` must be given in the command line of the executable. 
