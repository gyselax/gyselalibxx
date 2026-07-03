

# Class GMGPolarPoissonLikeSolver

**template &lt;class ToPhysicalMapping, class [**GridR**](structGridR.md), class [**GridTheta**](structGridTheta.md), class [**BSplinesR**](structBSplinesR.md), class [**BSplinesTheta**](structBSplinesTheta.md), class SplineBuilder, class SplineEvaluator&gt;**



[**ClassList**](annotated.md) **>** [**GMGPolarPoissonLikeSolver**](classGMGPolarPoissonLikeSolver.md)



_A Poisson-like solver using the GMGPolar multigrid library._ [More...](#detailed-description)

* `#include <gmg_polar_poisson_like_solver.hpp>`



Inherits the following classes: [IPolarPoissonLikeSolver](classIPolarPoissonLikeSolver.md)






















































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**GMGPolarPoissonLikeSolver**](#function-gmgpolarpoissonlikesolver) (ToPhysicalMapping to\_physical, SplineBuilder const & builder, SplineEvaluator const & evaluator, std::optional&lt; int &gt; max\_iterations=std::nullopt, std::optional&lt; double &gt; absTol=std::nullopt, std::optional&lt; double &gt; relTol=std::nullopt) <br>_Construct a_ [_**GMGPolarPoissonLikeSolver**_](classGMGPolarPoissonLikeSolver.md) _._ |
|  void | [**operator()**](#function-operator) (DField&lt; IdxRangeRTheta &gt; phi, DConstField&lt; IdxRangeRTheta &gt; rho) override const<br>_Solve the Poisson-like equation._  |
|  void | [**update\_coefficients**](#function-update_coefficients) (DConstField&lt; IdxRangeRTheta &gt; alpha, DConstField&lt; IdxRangeRTheta &gt; beta) override<br>_Rebuild the internal spline representations of α and β from grid values._  |
























































## Detailed Description


Solves -∇·(α∇φ) + βφ = ρ on a polar domain with homogeneous Dirichlet BCs at the outer boundary, using an across-the-origin discretisation at r = 0.




**Template parameters:**


* `ToPhysicalMapping` Mapping from (r,θ) to (x,y). 
* [**GridR**](structGridR.md) Discrete radial grid. 
* [**GridTheta**](structGridTheta.md) Discrete poloidal grid. 
* [**BSplinesR**](structBSplinesR.md) Radial B-spline space. 
* [**BSplinesTheta**](structBSplinesTheta.md) Poloidal B-spline space. 
* `SplineBuilder` 2D spline builder for ([**GridR**](structGridR.md) × [**GridTheta**](structGridTheta.md)). 
* `SplineEvaluator` 2D spline evaluator for ([**BSplinesR**](structBSplinesR.md) × [**BSplinesTheta**](structBSplinesTheta.md)). 




    
## Public Functions Documentation




### function GMGPolarPoissonLikeSolver 

_Construct a_ [_**GMGPolarPoissonLikeSolver**_](classGMGPolarPoissonLikeSolver.md) _._
```C++
inline GMGPolarPoissonLikeSolver::GMGPolarPoissonLikeSolver (
    ToPhysicalMapping to_physical,
    SplineBuilder const & builder,
    SplineEvaluator const & evaluator,
    std::optional< int > max_iterations=std::nullopt,
    std::optional< double > absTol=std::nullopt,
    std::optional< double > relTol=std::nullopt
) 
```





**Parameters:**


* `to_physical` The mapping from the logical to the physical domain. 
* `builder` A builder to construct the coefficients of the interpolation. 
* `evaluator` The evaluator for the interpolation. 
* `max_iterations` The maximum number of iterations that the solver should carry out. 
* `absTol` The absolute tolerance for the convergence of the solver. 
* `relTol` The relative tolerance for the convergence of the solver. 




        

<hr>



### function operator() 

_Solve the Poisson-like equation._ 
```C++
inline void GMGPolarPoissonLikeSolver::operator() (
    DField< IdxRangeRTheta > phi,
    DConstField< IdxRangeRTheta > rho
) override const
```





**Parameters:**


* `phi` The solution \(\phi\) on the grid. 
* `rho` The right-hand side \(\rho\) on the grid. 




        

<hr>



### function update\_coefficients 

_Rebuild the internal spline representations of α and β from grid values._ 
```C++
inline void GMGPolarPoissonLikeSolver::update_coefficients (
    DConstField< IdxRangeRTheta > alpha,
    DConstField< IdxRangeRTheta > beta
) override
```





**Parameters:**


* `alpha` Values of α at the grid interpolation points. 
* `beta` Values of β at the grid interpolation points. 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/pde_solvers/gmg_polar_poisson_like_solver.hpp`

