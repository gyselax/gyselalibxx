

# Class PolarSplineFEMPoissonLikeSolver

**template &lt;class [**GridR**](structGridR.md), class [**GridTheta**](structGridTheta.md), class [**PolarBSplinesRTheta**](structPolarBSplinesRTheta.md), concepts::Interpolation Interpolation2D, class Mapping, concepts::MetricTensor MetricTensorType, class IdxRangeFull&gt;**



[**ClassList**](annotated.md) **>** [**PolarSplineFEMPoissonLikeSolver**](classPolarSplineFEMPoissonLikeSolver.md)



_Define a polar PDE solver for a Poisson-like equation._ [More...](#detailed-description)

* `#include <polar_spline_fem_poisson_like_solver.hpp>`



Inherits the following classes: [IPolarPoissonLikeSolver](classIPolarPoissonLikeSolver.md)












## Classes

| Type | Name |
| ---: | :--- |
| class | [**CoeffEvaluator**](classPolarSplineFEMPoissonLikeSolver_1_1CoeffEvaluator.md) &lt;class Evaluator, class Coeff&gt;<br>_A wrapper that binds an evaluator with its batched coefficient field to present a single callable_ `double operator()(IdxBatch, CoordRTheta)` _._ |
| struct | [**QDimRMesh**](structPolarSplineFEMPoissonLikeSolver_1_1QDimRMesh.md) <br>_Tag the first dimension for the quadrature mesh._  |
| struct | [**QDimThetaMesh**](structPolarSplineFEMPoissonLikeSolver_1_1QDimThetaMesh.md) <br>_Tag the second dimension for the quadrature mesh._  |


## Public Types

| Type | Name |
| ---: | :--- |
| typedef typename GridR::continuous\_dimension\_type | [**R**](#typedef-r)  <br>_The radial dimension._  |
| typedef typename GridTheta::continuous\_dimension\_type | [**Theta**](#typedef-theta)  <br>_The poloidal dimension._  |








































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**PolarSplineFEMPoissonLikeSolver**](#function-polarsplinefempoissonlikesolver) (Mapping const & mapping, MetricTensorType const & metric\_tensor, Interpolation2D const & interpolation, IdxRangeFull idx\_range\_full, std::optional&lt; int &gt; max\_iter=std::nullopt, std::optional&lt; double &gt; res\_tol=std::nullopt, std::optional&lt; bool &gt; batch\_solver\_logger=std::nullopt, std::optional&lt; int &gt; preconditioner\_max\_block\_size=std::nullopt, std::optional&lt; GkoMatrixResidualNorm &gt; residual\_norm=std::nullopt) <br>_Instantiate a polar Poisson-like solver using FEM with B-splines._  |
|  void | [**operator()**](#function-operator) (PolarSplineRTheta spline, RHSFunction const & rhs) const<br>_Solve the Poisson-like equation._  |
|  void | [**operator()**](#function-operator_1) (DFieldBatchedRTheta phi, RHSFunction const & rhs) const<br>_Solve the Poisson-like equation._  |
|  void | [**operator()**](#function-operator_2) (DFieldBatchedRTheta phi, DConstFieldBatchedRTheta rhs) override const<br>_Solve the Poisson-like equation._  |
|  void | [**update\_coefficients**](#function-update_coefficients) (DConstFieldBatchedRTheta alpha, DConstFieldBatchedRTheta beta) override<br>_Update the coefficients_ \(alpha\) _and_\(beta\) _that define the equation._ |
























































## Detailed Description


Solve the following Partial Differential Equation


(1) \(L\phi = - \nabla \cdot (\alpha \nabla \phi) + \beta \phi = \rho\), in \(\Omega\),


\(\phi = 0\), on \(\partial \Omega\),


As finite element basis functions we will use polar b-splines which are divided into two types: 1) Basis splines that can be written as a tensor product of 1d basis splines ("non-singular B-splines") 2) Basis splines that cover the centre point and are defined as a linear combination of basis splines of type 1 ("singular B-splines")


(see in Emily Bourne's thesis "Non-Uniform Numerical Schemes for the Modelling of Turbulence
 in the 5D GYSELA Code". December 2022.)




**Template parameters:**


* [**GridR**](structGridR.md) The radial grid type. 
* [**GridTheta**](structGridTheta.md) The poloidal grid type. 
* [**PolarBSplinesRTheta**](structPolarBSplinesRTheta.md) The type of the 2D polar B-splines (on the coordinate system \((r,\theta)\) including B-splines which traverse the O point). 
* `Interpolation2D` The type of the 2D interpolation object used to evaluate the \(\alpha\) and \(\beta\) coefficient fields. Must satisfy concepts::Interpolation. 
* `Mapping` The type of the mapping from the logical domain to the physical domain where the equation is defined. 
* `MetricTensorType` The type of the metric tensor associated with the mapping. 
* `IdxRangeFull` The full index range of \(\phi\) including any batch dimensions. 




    
## Public Types Documentation




### typedef R 

_The radial dimension._ 
```C++
using PolarSplineFEMPoissonLikeSolver< GridR, GridTheta, PolarBSplinesRTheta, Interpolation2D, Mapping, MetricTensorType, IdxRangeFull >::R =  typename GridR::continuous_dimension_type;
```




<hr>



### typedef Theta 

_The poloidal dimension._ 
```C++
using PolarSplineFEMPoissonLikeSolver< GridR, GridTheta, PolarBSplinesRTheta, Interpolation2D, Mapping, MetricTensorType, IdxRangeFull >::Theta =  typename GridTheta::continuous_dimension_type;
```




<hr>
## Public Functions Documentation




### function PolarSplineFEMPoissonLikeSolver 

_Instantiate a polar Poisson-like solver using FEM with B-splines._ 
```C++
inline PolarSplineFEMPoissonLikeSolver::PolarSplineFEMPoissonLikeSolver (
    Mapping const & mapping,
    MetricTensorType const & metric_tensor,
    Interpolation2D const & interpolation,
    IdxRangeFull idx_range_full,
    std::optional< int > max_iter=std::nullopt,
    std::optional< double > res_tol=std::nullopt,
    std::optional< bool > batch_solver_logger=std::nullopt,
    std::optional< int > preconditioner_max_block_size=std::nullopt,
    std::optional< GkoMatrixResidualNorm > residual_norm=std::nullopt
) 
```



The equation we are studying:


(1) \(L\phi = - \nabla \cdot (\alpha \nabla \phi) + \beta \phi = \rho\), in \(\Omega\),


\(\phi = 0\), on \(\partial \Omega\).




**Parameters:**


* `mapping` The mapping from the logical domain to the physical domain where the equation is defined. 
* `metric_tensor` The metric tensor associated with the mapping. 
* `interpolation` An interpolation builder and evaluator to build coefficients allowing a function to be evaluated anywhere on the (r,theta) domain. 
* `idx_range_full` The full index range on which \(\phi\) is defined, including any batch dimensions ahead of [**GridR**](structGridR.md) and [**GridTheta**](structGridTheta.md). For the unbatched instantiation (IdxRangeFull == IdxRange&lt;GridR, GridTheta&gt;) this is simply the (r,theta) grid. 
* `max_iter` The maximum number of iterations possible for the batched CSR solver. 
* `res_tol` The residual tolerance for the batched CSR solver. The solver stops when \|\|r\|\|/\|\|b\|\| &lt;= res\_tol where r is the residual computed by the CG recurrence. 
* `batch_solver_logger` Indicates whether log information such as the residual and the number of iterations should be monitored. 
* `preconditioner_max_block_size` The maximum size of the Jacobi preconditioner used by the batched CSR solver. 
* `residual_norm` The norm used to evaluate the residual in the stopping criterion.



**Template parameters:**


* `Mapping` A class describing a mapping from curvilinear coordinates to Cartesian coordinates. 




        

<hr>



### function operator() 

_Solve the Poisson-like equation._ 
```C++
template<class RHSFunction>
inline void PolarSplineFEMPoissonLikeSolver::operator() (
    PolarSplineRTheta spline,
    RHSFunction const & rhs
) const
```



This operator returns the coefficients associated with the B-Splines of the solution \(\phi\).




**Parameters:**


* `spline` The spline representation of the solution \(\phi\), for every batch. 
* `rhs` The rhs \(\rho\) of the Poisson-like equation. The type is templated but we can use the [**PoissonLikeRHSFunction**](classPoissonLikeRHSFunction.md) class. It must be an object which can be called from GPU with an operator() which takes a batch index and a CoordRTheta. If there are no batch dimensions the operator() may take only a CoordRTheta. 




        

<hr>



### function operator() 

_Solve the Poisson-like equation._ 
```C++
template<class RHSFunction>
inline void PolarSplineFEMPoissonLikeSolver::operator() (
    DFieldBatchedRTheta phi,
    RHSFunction const & rhs
) const
```



This operator uses the other operator () and returns the values on the grid of the solution \(\phi\).




**Parameters:**


* `phi` The values of the solution \(\phi\) on the grid, for every batch. 
* `rhs` The rhs \(\rho\) of the Poisson-like equation. The type is templated but we can use the [**PoissonLikeRHSFunction**](classPoissonLikeRHSFunction.md) class. It must be an object which can be called from GPU with an operator() which takes a batch index and a CoordRTheta. If there are no batch dimensions the operator() may take only a CoordRTheta. 




        

<hr>



### function operator() 

_Solve the Poisson-like equation._ 
```C++
inline void PolarSplineFEMPoissonLikeSolver::operator() (
    DFieldBatchedRTheta phi,
    DConstFieldBatchedRTheta rhs
) override const
```



This operator uses the other operator () and returns the values on the grid of the solution \(\phi\).




**Parameters:**


* `phi` The values of the solution \(\phi\) on the grid, for every batch. 
* `rhs` The rhs \(\rho\) of the Poisson-like equation on the grid, for every batch. 




        

<hr>



### function update\_coefficients 

_Update the coefficients_ \(alpha\) _and_\(beta\) _that define the equation._
```C++
inline void PolarSplineFEMPoissonLikeSolver::update_coefficients (
    DConstFieldBatchedRTheta alpha,
    DConstFieldBatchedRTheta beta
) override
```





**Parameters:**


* `alpha` The \(\alpha\) function in the definition of the Poisson-like equation defined at the grid points. 
* `beta` The \(\beta\) function in the definition of the Poisson-like equation defined at the grid points. 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/pde_solvers/polar_spline_fem_poisson_like_solver.hpp`

