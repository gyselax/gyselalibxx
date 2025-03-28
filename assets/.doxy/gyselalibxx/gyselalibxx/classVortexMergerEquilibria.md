

# Class VortexMergerEquilibria

**template &lt;class Mapping&gt;**



[**ClassList**](annotated.md) **>** [**VortexMergerEquilibria**](classVortexMergerEquilibria.md)



_Equilibrium solution of a Vlasov-Poissson equations system in polar coordinates._ [More...](#detailed-description)

* `#include <vortex_merger_equilibrium.hpp>`





































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**VortexMergerEquilibria**](#function-vortexmergerequilibria) (Mapping const & mapping, IdxRangeRTheta const & grid, SplineRThetaBuilder\_host const & builder, SplineRThetaEvaluatorNullBound\_host const & evaluator, [**PolarSplineFEMPoissonLikeSolver**](classPolarSplineFEMPoissonLikeSolver.md)&lt; [**GridR**](structGridR.md), [**GridTheta**](structGridTheta.md), [**PolarBSplinesRTheta**](structPolarBSplinesRTheta.md), SplineRThetaEvaluatorNullBound &gt; const & poisson\_solver) <br>_Instantiate a_ [_**VortexMergerEquilibria**_](classVortexMergerEquilibria.md) _._ |
|  void | [**find\_equilibrium**](#function-find_equilibrium) (host\_t&lt; DFieldRTheta &gt; sigma, host\_t&lt; DFieldRTheta &gt; phi\_eq, host\_t&lt; DFieldRTheta &gt; rho\_eq, std::function&lt; double(double const)&gt; const & function, double const phi\_max, double const tau, int count\_max=25) const<br>_Get an equilibrium._  |
|  void | [**set\_equilibrium**](#function-set_equilibrium) (host\_t&lt; DFieldRTheta &gt; rho\_eq, std::function&lt; double(double const)&gt; function, double const phi\_max, double const tau) <br>_Set an equilibrium._  |




























## Detailed Description




**Template parameters:**


* `Mapping` A class describing a mapping from curvilinear coordinates to Cartesian coordinates. 




    
## Public Functions Documentation




### function VortexMergerEquilibria 

_Instantiate a_ [_**VortexMergerEquilibria**_](classVortexMergerEquilibria.md) _._
```C++
inline VortexMergerEquilibria::VortexMergerEquilibria (
    Mapping const & mapping,
    IdxRangeRTheta const & grid,
    SplineRThetaBuilder_host const & builder,
    SplineRThetaEvaluatorNullBound_host const & evaluator,
    PolarSplineFEMPoissonLikeSolver < GridR , GridTheta , PolarBSplinesRTheta , SplineRThetaEvaluatorNullBound > const & poisson_solver
) 
```





**Parameters:**


* `mapping` The mapping function from the logical index range to the physical index range. 
* `grid` The index range where the equilibrium is defined. 
* `builder` A spline builder to get the spline representation of the RHS of the PDE. 
* `evaluator` The evaluator of B-splines for the RHS of the PDE. 
* `poisson_solver` The PDE solver which computes the electrical potential. 




        

<hr>



### function find\_equilibrium 

_Get an equilibrium._ 
```C++
inline void VortexMergerEquilibria::find_equilibrium (
    host_t< DFieldRTheta > sigma,
    host_t< DFieldRTheta > phi_eq,
    host_t< DFieldRTheta > rho_eq,
    std::function< double(double const)> const & function,
    double const phi_max,
    double const tau,
    int count_max=25
) const
```



The equilibrium is determined by the eigenvalue problem. For the given initial data ,



* compute ;
* compute  with ;
* for  given, compute  with  ;
* compute .




We iterate until .


For the vortex merger simulation, .


The algorithm is also detailed in Edoardo Zoni's article ([https://doi.org/10.1016/j.jcp.2019.108889](https://doi.org/10.1016/j.jcp.2019.108889)).




**Parameters:**


* `sigma` The  parameter. 
* `phi_eq` The equilibrium electrical potential . 
* `rho_eq` The equilibrium density . 
* `function` The function . 
* `phi_max` The maximal value of the electrical potential . 
* `tau` The  parameter. 
* `count_max` The maximal number of iteration of the implicit loop. 




        

<hr>



### function set\_equilibrium 

_Set an equilibrium._ 
```C++
inline void VortexMergerEquilibria::set_equilibrium (
    host_t< DFieldRTheta > rho_eq,
    std::function< double(double const)> function,
    double const phi_max,
    double const tau
) 
```





**Parameters:**


* `rho_eq` The equilibrium density . 
 
* `function` The function  used to compute the equilibrium. 
* `phi_max` The maximal value of the electrical potential . 
* `tau` The  parameter. 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/geometryRTheta/initialisation/vortex_merger_equilibrium.hpp`

