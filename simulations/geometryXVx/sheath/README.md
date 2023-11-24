# Sheath simulations

## Description 
The `sheath` executable allows the study of kinetic plasma-wall interaction. A Boltzmann-Poisson system is solved for the electric field and the distribution function for both electrons and ions species. The model is one dimensional in space and velocity (1D1V).   

The particularity of the model lies in the description of the wall region, which is treated using immersed boundary condition (or penalization). Electrons and ions are absorbed inside this region. The inertia difference between the two species drives the appearance of a positively charged layer (the sheath) at the plasma boundary. More details about the sheath physics along with the need for non-equidistant meshes to conduct such simulations can be found in [1,2].  

## Usage
After building the code, run the executable located in `build/simulations/sheath/`. To use the default simulation parameters the user can provide the `--dump-config` option when launching the executable.

## Recommended parameters
Two sets of parameters are available : 
- the default parameter given in `sheath.yaml.hpp`. This cas corresponds to a very light simulation case that runs a few iterations for testing purposes. This is the set of parameters that can be retrieved with the `--dump-config` command. 
- the parameters given in the folder `ref_simulation`, in the `sheath_ref.yaml` file. This simulation presents a glimpse of the plasma-wall interaction physics that can be seen on the figures in the mentioned folder: formation of a positively charged layer in front of the wall region, accompanied by a drop of the electric potential that confines fast electrons. A particle flux going towards the wall exists in the plasma. With the number of iterations performed to obtain these figures, the system is not at steady state. To observe a steady state plasma-wall interaction with interesting physical features (supersonic flow of ions, truncation of fast velocity electrons, etc.) see the paramerers of [2].

## Verification of the simulation
The accuracy of the physical results from a sheath simulation can be verified as follows. The simulation should be run with the parameters given in `sheath.yaml.hpp` file. Then, any conservation diagnostic can be used to test the accuracy of the results. For instance, a post-process script that breaks down the terms present in the equation of energy conservation is `post-process/PythonScripts/geometryXVx/sheath/plot_conservation_energy`. This script can be run as an executable in the folder containing the simulation results. In the simple XVx geometry The equation of energy conservation is written as 

$\partial_t \Pi_a + 1/\sqrt{A_a}\, \partial_x Q_a - 2\Gamma_a q_a E/\sqrt{A_a} = S_{e}$  

Where $\Gamma_a$, $\Pi_a$ and $Q_a$ stand for the particle flux, momentum flux and energy flux of species $a$ (typically ions or electrons) respectively. These quantity are computed by the script using the distribution function of each species which is an output of the code. $A_a = m_e/m_a$ is the mass ratio of species $a$, and $q_a$ represents its charge ($+e$ for ions, $-e$ for electrons with $e$ the elementary charge). $E$ is the electric field, and $S_{e}$ the energy source term. The script plots all of the terms that appear in the energy conservation equation, along with the error (denoted as `lhs-rhs`) defined as being the left-hand-side term minus the right-hand-side terms, i.e.
$lhs-rhs = \partial_t \Pi_a + 1/\sqrt{A_a}\, \partial_x Q_a - 2\Gamma_a q_a E/\sqrt{A_a} - S_{e}$. This error term is plotted as well, and thus a simulation can be verified by ensuring that this error remains small as compared to the other terms. It is recommended to look at the conservation equation for electrons first, since as the lightest species, the numerical error is larger. The script also produces a 2D map of the error as a function of time and space. One shall note that in order to compute the time derivative term $\partial_t \Pi_a$, it is important to save all of the timesteps of the simulations, i.e. the parameters `time_diag` and `deltat` should be equal, as in the `sheath.yaml.hpp` file. Lastly, if one wishes to test a simulation that is faster, the number of iterations can be reduced to 10 for instance. The folder `conservation_plots` contains the graphs of the conservation equations computed using a simulation with the parameters of the `sheath.yaml.hpp` file.

## References
- [1] E. Bourne, Y. Munschy, V. Grandgirard, M. Mehrenberger, and P. Ghendrih, Non-Uniform Splines for Semi-Lagrangian Kinetic Simulations of the Plasma Sheath (2022)
- [2] Y. Munschy, E. Bourne, P. Ghendrih, G. Dif-Pradalier, Y. Sarazin, V. Grandgirard, and P. Donnel, Kinetic plasma-wall interaction using immersed boundary conditions (2023)
