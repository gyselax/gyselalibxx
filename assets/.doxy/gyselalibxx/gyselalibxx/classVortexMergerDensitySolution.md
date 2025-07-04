

# Class VortexMergerDensitySolution

**template &lt;class Mapping&gt;**



[**ClassList**](annotated.md) **>** [**VortexMergerDensitySolution**](classVortexMergerDensitySolution.md)



_Initial condition for the vortex merger simulation._ [More...](#detailed-description)

* `#include <vortex_merger_initialisation.hpp>`





































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**VortexMergerDensitySolution**](#function-vortexmergerdensitysolution) (Mapping const & mapping) <br>_Instantiate a_ [_**VortexMergerDensitySolution**_](classVortexMergerDensitySolution.md) _._ |
|  void | [**set\_initialisation**](#function-set_initialisation) (host\_t&lt; DFieldRTheta &gt; rho\_init, host\_t&lt; DConstFieldRTheta &gt; rho\_eq, const double eps, const double sigma, const double x\_star\_1, const double y\_star\_1, const double x\_star\_2, const double y\_star\_2) <br>_Set an initial condition for the vortex merger simulation._  |




























## Detailed Description




**Template parameters:**


* `Mapping` A class describing a mapping from curvilinear coordinates to Cartesian coordinates. 




    
## Public Functions Documentation




### function VortexMergerDensitySolution 

_Instantiate a_ [_**VortexMergerDensitySolution**_](classVortexMergerDensitySolution.md) _._
```C++
inline explicit VortexMergerDensitySolution::VortexMergerDensitySolution (
    Mapping const & mapping
) 
```





**Parameters:**


* `mapping` A mapping function from the logical domain to the physical domain. 




        

<hr>



### function set\_initialisation 

_Set an initial condition for the vortex merger simulation._ 
```C++
inline void VortexMergerDensitySolution::set_initialisation (
    host_t< DFieldRTheta > rho_init,
    host_t< DConstFieldRTheta > rho_eq,
    const double eps,
    const double sigma,
    const double x_star_1,
    const double y_star_1,
    const double x_star_2,
    const double y_star_2
) 
```



The initial condition is given by  
\[\rho(0, x, y) = \rho_{eq}(x,y) + \varepsilon
     \left(
         \exp\left[-\frac{(x - x_1^*)^2 + (y - y_1^*)^2}{2\sigma^2}\right]
      +  \exp\left[-\frac{(x - x_2^*)^2 + (y - y_2^*)^2}{2\sigma^2}\right]
     \right)\]



with \(\rho_{eq}\) given by [**VortexMergerEquilibria::set\_equilibrium**](classVortexMergerEquilibria.md#function-set_equilibrium).


The initial condition is also given in Edoardo Zoni's article ([https://doi.org/10.1016/j.jcp.2019.108889](https://doi.org/10.1016/j.jcp.2019.108889)).




**Parameters:**


* `rho_init` The initial condition of the density. 
* `rho_eq` The equilibrium density solution computed thanks to [**VortexMergerEquilibria::set\_equilibrium**](classVortexMergerEquilibria.md#function-set_equilibrium). 
* `eps` The \(\varepsilon\) amplitude of the perturbation. 
* `sigma` The \(\sigma\) of the Gaussian functions. 
* `x_star_1` The \(x_1^*\) parameter. 
* `y_star_1` The \(y_1^*\) parameter. 
* `x_star_2` The \(x_2^*\) parameter. 
* `y_star_2` The \(y_2^*\) parameter. 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/geometryRTheta/initialisation/vortex_merger_initialisation.hpp`

