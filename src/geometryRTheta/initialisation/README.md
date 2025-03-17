# Initialisation

This folder defines initialisation and equilibrium functions for the simulations.  
For more details on each simulation, see [simulations](./../../../simulations/README.md).

## Diocotron instability

More details about the diocotron simulation are given in [diocotron](./../../../simulations/geometryRTheta/diocotron/README.md).

The initializatrion and equilibrium functions are defined in the DiocotronDensitySolution class. 

### Initialisation

We suppose as perturbed initial condition
 - \(\rho(0,r, \theta) = ( 1 + \epsilon \cos(l\theta)) \ \mathbb{I}_{ R_1 \leq r \leq R_2 }\)

where \(\mathbb{I}\) is the characteristic function.

In Edoardo Zoni's article [1], an exponential is also added to make the solution smoother: 
- \(\rho(0,r, \theta) = ( 1 + \epsilon \cos(l\theta)) \exp\left(- \left(\frac{r - \bar{r}}{d}\right)^p\right) \ \mathbb{I}_{ R_1 \leq r \leq R_2 }\)  

with $p = 50$, \(\bar{r} = \frac{R_1 + R_2}{2}\) and \(d = \frac{R_2 - R_1}{2}\). 
This last version is implemented in the code. 

### Equilibrium
The associated equilibrium is 
 - \(\rho_0(r) = \mathbb{I}_{ R_1 \leq r \leq R_2 }\)

or 
 - \(\rho_0(r) = \exp\left(- \left(\frac{r - \bar{r}}{d}\right)^p\right) \ \mathbb{I}_{ R_1 \leq r \leq R_2 }\)
in  Edoardo Zoni's article [1].


### Dispersion relation

The DiocotronDensitySolution also solves the dispersion relation: 

```math
\left(\frac{\omega}{\omega_D} \right)^2 - b_l \frac{\omega}{\omega_D} + c_l = 0
```

with 
```math
b_l \left( 1 - \frac{W_1^{2l}}{W_2^{2l}} \right) = l\left[ 1 - \frac{R_1^2}{R_2^2} + \frac{\omega_q }{\omega_D} \left( 1 + \frac{R_1^2}{R_2^2} \right)\right]\left[ 1 - \frac{W_1^{2l}}{W_2^{2l}} \right] + \left(1 - \frac{R_1^{2l}}{R_2^{2l}} \right) \left(\frac{R_2^{2l}}{W_2^{2l}} - \frac{W_1^{2l}}{R_1^{2l}} \right)
```

```math
 c_l \left( 1 - \frac{W_1^{2l}}{W_2^{2l}} \right) = l^2 \frac{\omega_q }{\omega_D} \left[1 - \frac{R_1^2}{R_2^2} + \frac{\omega_q }{\omega_D} \frac{R_1^2}{R_2^2}\right]\left( 1 - \frac{W_1^{2l}}{W_2^{2l}}  \right) -l \frac{\omega_q }{\omega_D} \left( 1 - \frac{W_1^{2l}}{R_2^{2l}} \right)\left( 1 - \frac{R_2^{2l}}{W_2^{2l}} \right) 
\\ \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad + l \left[1 - \frac{R_1^2}{R_2^2} + \frac{\omega_q }{\omega_D} \frac{R_1^2}{R_2^2} \right]\left( 1 - \frac{R_1^{2l}}{W_2^{2l}} \right) \left( 1 - \frac{W_1^{2l}}{R_1^{2l}} \right) - \left( 1 - \frac{R_2^{2l}}{W_2^{2l}} \right) \left( 1 - \frac{W_1^{2l}}{R_1^{2l}} \right)  \left( 1 - \frac{R_1^{2l}}{R_2^{2l}} \right)
```
 
The imaginary part of the dispersion relation solution is given by DiocotronDensitySolution::get\_slope() 
and its real part by DiocotronDensitySolution::get\_frequency(). The imaginary part of the dispersion relation solution 
corresponds to the slope of the L2-norm perturbation in the linear phase ($t \in [20, 50]$ s for the default parameters).
 


## Vortex merger

More details about the diocotron simulation are given in [vortex\_merger](./../../../simulations/geometryRTheta/vortex_merger/README.md).
 
### Equilibrium
The equilibrium is determined by the eigenvalue problem of finding $(\sigma, \phi)$ such that 

```math
    - \nabla \cdot \nabla \phi = \sigma f(\phi),
```

with given $f$  such that $f'(\phi) \neq 0$. See Edoardo Zoni's article [1] for more details. 

The algorithm for given initial data $(\sigma^0, \phi^0)$ is the following: 
* compute the density: $\rho^{i} = \sigma^{i-1} f(\phi^{i-1})$, 
* compute a temporary electrical potential: \(- \nabla\cdot\nabla \phi_*^{i} = \rho^i\), 
* compute the coefficient: 
    * if the maximum value \(\phi_{\text{max}}\) is given: \(c^i = \frac{\phi_{\text{max}}}{\Vert\phi_*^i\Vert_{\mathcal{L}^\infty}}\), 
* update the data: \((\sigma^i, \phi^i) = c^i (\sigma^{i-1}, \phi_*^i)\). 

Then, we repeat all the steps until $|\sigma^i - \sigma^{i-1}| \leq \tau$ for a given tolerance. 

For the vortex merger simulation, we choose \(f(\phi)= \phi^2\) and \(\phi_{\text{max}} = 1\). 


### References
[1]    Edoardo Zoni, Yaman Güçlü, "Solving hyperbolic-elliptic problems on singular mapped disk-like domains with the 
method of characteristics and spline finite elements", https://doi.org/10.1016/j.jcp.2019.108889, Journal of Computational Physics, 2019.

## Contents 

* diocotron\_initialisation\_equilibrium.hpp : defines a DiocotronDensitySolution class which computes the initialisation and equilibrium functions. 
* vortex\_merger\_equilibrium.hpp : defines a VortexMergerEquilibria class which computes the equilibrium function. 
* vortex\_merger\_initialisation.hpp : defines a VortexMergerDensitySolution class which computes the initialisation function. 


