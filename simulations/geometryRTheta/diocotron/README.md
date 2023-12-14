# Diocotron instability

## Studied problem
### Field setup 
We work with cylindrical coordinates $`(r,\theta,z)`$ on the corresponding orthonormal basis vectors $`(e_r, e_\theta, e_z)`$.
We suppose that the particles (electrons) are trapped between two cylindrically conducting walls at the radii $`W_1`$ and $`W_2`$, (setup from Petri's articles [6], [7]). 
The plasma is confined between the radii $`R_1`$ et $`R_2`$ such that 

```math
W_1 \leq R_1 \lt R_2 \leq W_2. 
```


### General case
Conservation of the density equation and Maxwell-Poisson equation:

```math
\\ \partial_t \rho + div(\rho v) = 0 
\\ \Delta \phi + q\frac{\rho}{\varepsilon_0} = 0 
\\ v = \frac{E\wedge B}{|B|^2} 
\\ E = -\nabla \phi 
```

where
* $`\rho = \rho(t,r,\theta)`$ is the **density** of the particles (electrons with mass $`m_e`$ and charge $q = -1$);
* $v$ is the velocity of the particles ($div(v) = 0$); 
* $\phi$ is the usual electric potential associated withe the total electric field; 
* $(E,B)$ the electromagnetic field such that $`B = B_z(r) e_z`$ the external magnetic field


### Simplified case
We study a simplified case where 
- $`B = e_z`$, 
- $`\varepsilon_0 = 1`$. 

On a circular mapping, it gives 
$`v = -\nabla\phi \wedge e_z`$, $div(v) = 0$  
and 
$`v = -\nabla\phi \wedge e_z = - \frac{\partial_\theta\phi}{r^2}e_r + \partial_r\phi e_\theta`$ with $`(e_r, e_\theta)`$ 
the unnormalized local contravariant base, 

```math
 \partial_t \rho - \frac{\partial_\theta\phi}{r}\partial_r\rho + \partial_r\phi \frac{\partial_\theta\rho}{r} =0, 
```

```math
 \Delta \phi = \rho. 
```


## Test case - diocotron instability
The test case implemented tests the diocoton instabilities (see Davidson book's [3], Edoardo Zoni's article [2] and article [4]).
The initial function and equilibrium are defined in DiocotronDensitySolution and documented in 
[initialization](./../../../src/geometryRTheta/initialization/README.md). 

We suppose as perturbed initial condition

```math
\rho(0,r, \theta) = ( 1 + \epsilon \cos(l\theta)) \ \mathbb{I}_{ R_1 \leq r \leq R_2 }
```

where $\mathbb{I}$ is the characteristic function.

In Edoardo Zoni's article [2], an exponential is also added to make the solution smoother: 

```math
\rho(0,r, \theta) = ( 1 + \epsilon \cos(l\theta)) \exp\left(- \left(\frac{r - \bar{r}}{d}\right)^p\right) \ \mathbb{I}_{ R_1 \leq r \leq R_2 }
```

with $`p = 50 \text{, } \bar{r} = \frac{R_1 + R_2}{2} \text{, and } d = \frac{R_2 - R_1}{2}`$. 

The explicit solution $(\phi,\rho)$ is given by 
```math
 \rho(t,r, \theta) = \rho_{0}(r) + \varepsilon\rho_{1}(t,r,\theta),
 \phi(t,r, \theta) = \phi_{0}(r) + \varepsilon\phi_{1}(t,r,\theta),
```

with 
```math
\rho_0(r) = \mathbb{I}_{ R_1 \leq r \leq R_2 }
```

or 

```math
\rho_0(r) = \exp\left(- \left(\frac{r - \bar{r}}{d}\right)^p\right) \ \mathbb{I}_{R_1 \leq r \leq R_2}
```

 in  Edoardo Zoni's article [2], and
 ```math
 \rho_{1}(t,r, \theta) = \widehat{\rho_{1,l}}(r) e^{il\theta} e^{-i\omega t},
 \phi_{1}(t,r, \theta) = \widehat{\phi_{1,l}}(r) e^{il\theta} e^{-i\omega t},
```



with $`l\in\mathbb{N}`$ the selected mode, $`\omega \in \mathbb{C}`$ a pulsation to be determined, and 
 - if  $`W_1 \leq r\lt R_1`$ then,
   $`\widehat{\phi_{1,l}} (r) = \phi_{1,I}(r) = \frac{R_1^l}{r^l} \left( B R_1^l + C R_1^{-l}\right) \frac{r^{2l} - W_1^{2l}}{R_1^{2l} - W_1^{2l}}`$
 - if  $`R_1 \leq r\lt R_2`$ then,
   $`\widehat{\phi_{1,l}} (r) = \phi_{1,II}(r) =  B r^l + C r^{-l}`$ 
 - if  $`R_2 \leq r\lt W_2`$ then,
   $`\widehat{\phi_{1,l}} (r) = \phi_{1,III}(r) = \frac{R_2^l}{r^l} \left( B R_2^l + C R_2^{-l}\right) \frac{r^{2l} - W_2^{2l}}{ R_2^{2l} - W_2^{2l}}`$.

for $B$ and $C$ constants stastifying continuity conditions such that $`A[B,C]^t = 0`$ 

with 
```math
 A_{1,1} = R_1^l -  \frac{R_1^{2l+1}}{R_1^{2l} - W_1^{2l}} \left[  R_1^{l-1} + \frac{1}{R_1^{l+1}}  W_1^{2l}  \right] + \frac{ R_1^l}{\omega - l\omega_m(R_1)},
\\ A_{1,2} = - R_1^{-l} - \frac{R_1}{R_1^{2l} - W_1^{2l}} \left[  R_1^{l-1} + \frac{1}{R_1^{l+1}}  W_1^{2l}  \right] +  \frac{ R_1^{-l}}{\omega - l\omega_m(R_1)},
\\ A_{2,1} = - R_2^{l} + \frac{R_2^{2l+1}}{R_2^{2l} - W_2^{2l}} \left[  R_2^{l-1} + \frac{1}{R_2^{l+1}}  W_2^{2l}  \right] - \frac{ R_2^l }{\omega - l\omega_m(R_2)},
\\ A_{2,2} = R_2^{-l} + \frac{R_2}{R_2^{2l} - W_2^{2l}} \left[  R_2^{l-1} + \frac{1}{R_2^{l+1}}  W_2^{2l}  \right] - \frac{ R_2^{-l} }{\omega - l\omega_m(R_2)}.
```

where $\omega$ satifies 
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

and $`\omega_q = - \frac{2Q}{R_1^2}`$, $`\omega_D = \frac{1}{2}`$ and $Q$ the charge carried by unit of length at $`r = W_1`$ 
(in our simulation, we put $`Q = 0`$).
 
 
 
 ### Predictor-corrector
 
 To solve the equations system, we use a predictor-corrector methods. Severals predictor-corrector methods are
 implemented in ITimeSolverRP and documented in [time\_solver](./../../../src/geometryRTheta/time_solver/README.md).


## References
[1]    Emily Bourne. “Non-Uniform Numerical Schemes for the Modelling of Turbulence in the 5D GY-SELA Code”. PhD thesis. Aix-Marseille Université, Dec. 2022.

[2]    Edoardo Zoni, Yaman Güçlü, "Solving hyperbolic-elliptic problems on singular mapped disk-like domains with the 
method of characteristics and spline finite elements", https://doi.org/10.1016/j.jcp.2019.108889, Journal of Computational Physics, 2019.

[3]    R. C. Davidson.Physics of non neutral plasmas. 1990. Chap. 6, The Dioctron instability.

[4]    Eric Madaule, Sever Adrian Hirstoaga, Michel Mehrenberger and Jerôme Petri. 
“Semi-Lagrangian simulations of the diocotron instability”. HAL (July 2013). doi: https://hal.inria.fr/hal-00841504.

[5]    Eric Sonnendrucker, Jean Roche, Pierre Bertrand and Alain Gizzo. 
“The Semi-Lagrangian Method for the Numerical Resolution of the Vlasov Equation”. Journal of Computational Physics (1999).

[6]    Jerôme  Petri.  “Non-linear  evolution  of  the  diocotron  instability  in  a  pulsar  electrosphere: two-dimensional particle-in-cell simulations”. Astronomy & Astrophysics (Apr. 2009).

[7]    Jerôme  Petri. “The diocotron instability in a pulsar ”cylindrical” electrosphere”. Astronomy & Astrophysics(Nov. 2006)


## Contents

 - diocotron.cpp : runs a diocotron instability simulation.
 - params.yaml : defines the parameters of the simulation. 
 
 Python files: (in `post_process/PythonScripts/geometryRTheta/diocotron` folder)
 - animation\_rho\_phi.py : plot and save the density and electrical potential function in time. 
 - display\_L2\_norms.py : plot the L2 norms of the perturbation of the density and electrical potential function in time. 
 - mass\_conservation.py : plot the relative errors of the mass of the particles. 
 
 
 ### Recommended parameters: 
* `Mesh:`
  * `r_size: 128` : number of cells in $r$-dimension. (Tests in Edoardo Zoni's article.)
  * `p_size: 256` : number of cells in $\theta$-dimension. (Tests in Edoardo Zoni's article.)
  * `r_min: 0.0`  : position of the inner wall $`W_1`$. (Tests in Edoardo Zoni's article.)
  * `r_minus: 0.45` : position of the inner boundary of the initial density $`R_1`$. (Tests in Edoardo Zoni's article.)
  * `r_plus: 0.50`: position of the outer boundary of the initial density $`R_2`$. (Tests in Edoardo Zoni's article.)
  * `r_max: 1.0` : position of the outer wall $`W_2`$. (Tests in Edoardo Zoni's article.)

* `Time:`
  * `delta_t: 0.1` : time step. (Tests in Edoardo Zoni's article.)
  * `final_T: 70.0`: final time of the simulation (end of the linear phase at $T = 50s$). (Tests in Edoardo Zoni's article.)
  
* `Perturbation:`
  * `charge_Q: 0` : charge carried by the inner conductor at $`r = W_1`$.
  * `l_mode: 9` : mode of the pertubation $\varepsilon \cos(lx)$. 
  * `eps: 0.0001` : amplitude of the pertubation $\varepsilon \cos(lx)$. 
  
* `Output:`
  * `time_step_diag: 10` : number of time steps between two recordings of the data. 
 
 
 ### Executables:
* `diocotron_EXPLICIT_PREDCORRR_EULER_METHOD` : uses the predictor-corrector defined in BslExplicitPredCorrRP (ill-defined for other time integration methods). 
* `diocotron_IMPLICIT_PREDCORRR_EULER_METHOD` : uses the predictor-corrector defined in BslImplicitPredCorrRP (ill-defined for other time integration methods).
* `diocotron_PREDCORRR_EULER_METHOD` : uses the predictor-corrector defined in BslPredCorrRP with a Euler method for the BslAdvectionRP advection operator.
* `diocotron_PREDCORRR_CRANK_NICOLSON_METHOD` : uses the predictor-corrector defined in BslPredCorrRP with a CrankNicolson method for the BslAdvectionRP advection operator.
* `diocotron_PREDCORRR_RK3_METHOD` : uses the predictor-corrector defined in BslPredCorrRP with a RK3 method for the BslAdvectionRP advection operator.
* `diocotron_PREDCORRR_RK4_METHOD` : uses the predictor-corrector defined in BslPredCorrRP with a RK4 method for the BslAdvectionRP advection operator.
 
The results of the simulation are saved in a `output` folder that the executable creates. 

The path to the `params.yaml` must be given in the command line of the executable. 




 
