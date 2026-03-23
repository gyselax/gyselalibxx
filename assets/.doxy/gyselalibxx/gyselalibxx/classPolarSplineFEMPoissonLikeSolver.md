

# Class PolarSplineFEMPoissonLikeSolver

**template &lt;class [**GridR**](structGridR.md), class [**GridTheta**](structGridTheta.md), class [**PolarBSplinesRTheta**](structPolarBSplinesRTheta.md), class SplineRThetaEvaluatorNullBound, class Mapping, class IdxRangeFull&gt;**



[**ClassList**](annotated.md) **>** [**PolarSplineFEMPoissonLikeSolver**](classPolarSplineFEMPoissonLikeSolver.md)



_Define a polar PDE solver for a Poisson-like equation._ [More...](#detailed-description)

* `#include <polarpoissonlikesolver.hpp>`















## Classes

| Type | Name |
| ---: | :--- |
| struct | [**InternalBatchDim**](structPolarSplineFEMPoissonLikeSolver_1_1InternalBatchDim.md) <br>_The tag for the batch dimension for the equation. This is public due to Cuda._  |
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
|   | [**PolarSplineFEMPoissonLikeSolver**](#function-polarsplinefempoissonlikesolver) (Mapping const & mapping, SplineRThetaEvaluatorNullBound const & spline\_evaluator, std::optional&lt; int &gt; max\_iter=std::nullopt, std::optional&lt; double &gt; res\_tol=std::nullopt, std::optional&lt; bool &gt; batch\_solver\_logger=std::nullopt, std::optional&lt; int &gt; preconditioner\_max\_block\_size=std::nullopt) <br>_Instantiate a polar Poisson-like solver using FEM with B-splines._  |
|  void | [**operator()**](#function-operator) (RHSFunction const & rhs, PolarSplineRTheta spline) const<br>_Solve the Poisson-like equation._  |
|  void | [**operator()**](#function-operator_1) (RHSFunction const & rhs, DFieldRTheta phi) const<br>_Solve the Poisson-like equation._  |
|  void | [**update\_coefficients**](#function-update_coefficients) (ConstSpline2D coeff\_alpha, ConstSpline2D coeff\_beta) <br>_Update the coefficients_ \(alpha\) _and_\(beta\) _that define the equation._ |


## Public Static Functions

| Type | Name |
| ---: | :--- |
|  KOKKOS\_INLINE\_FUNCTION double | [**get\_polar\_bspline\_vals**](#function-get_polar_bspline_vals) (CoordRTheta coord, IdxBSPolar idx) <br>_Get the value of the specified polar bspline at the specified point._  |


























## Detailed Description


Solve the following Partial Differential Equation


(1) \(L\phi = - \nabla \cdot (\alpha \nabla \phi) + \beta \phi = \rho\), in \(\Omega\),


\(\phi = 0\), on \(\partial \Omega\),


As finite element basis functions we will use polar b-splines which are divided into two types: 1) Basis splines that can be written as a tensor product of 1d basis splines ("non-singular B-splines") 2) Basis splines that cover the centre point and are defined as a linear combination of basis splines of type 1 ("singular B-splines")


(see in Emily Bourne's thesis "Non-Uniform Numerical Schemes for the Modelling of Turbulence
 in the 5D GYSELA Code". December 2022.)




**Template parameters:**


* [**GridR**](structGridR.md) The radial grid type. 
* [**GridR**](structGridR.md) The poloidal grid type. 
* [**PolarBSplinesRTheta**](structPolarBSplinesRTheta.md) The type of the 2D polar B-splines (on the coordinate system \((r,\theta)\) including B-splines which traverse the O point). 
* `SplineRThetaEvaluatorNullBound` The type of the 2D (cross-product) spline evaluator. 
* `Mapping` The type of the mapping from the logical domain to the physical domain where the equation is defined. 
* `IdxRangeFull` The full index range of \(\phi\) including any batch dimensions. 




    
## Public Types Documentation




### typedef R 

_The radial dimension._ 
```C++
using PolarSplineFEMPoissonLikeSolver< GridR, GridTheta, PolarBSplinesRTheta, SplineRThetaEvaluatorNullBound, Mapping, IdxRangeFull >::R =  typename GridR::continuous_dimension_type;
```




<hr>



### typedef Theta 

_The poloidal dimension._ 
```C++
using PolarSplineFEMPoissonLikeSolver< GridR, GridTheta, PolarBSplinesRTheta, SplineRThetaEvaluatorNullBound, Mapping, IdxRangeFull >::Theta =  typename GridTheta::continuous_dimension_type;
```




<hr>
## Public Functions Documentation




### function PolarSplineFEMPoissonLikeSolver 

_Instantiate a polar Poisson-like solver using FEM with B-splines._ 
```C++
inline PolarSplineFEMPoissonLikeSolver::PolarSplineFEMPoissonLikeSolver (
    Mapping const & mapping,
    SplineRThetaEvaluatorNullBound const & spline_evaluator,
    std::optional< int > max_iter=std::nullopt,
    std::optional< double > res_tol=std::nullopt,
    std::optional< bool > batch_solver_logger=std::nullopt,
    std::optional< int > preconditioner_max_block_size=std::nullopt
) 
```



The equation we are studying:


(1) \(L\phi = - \nabla \cdot (\alpha \nabla \phi) + \beta \phi = \rho\), in \(\Omega\),


\(\phi = 0\), on \(\partial \Omega\).




**Parameters:**


* `mapping` The mapping from the logical domain to the physical domain where the equation is defined. 
* `spline_evaluator` An evaluator for evaluating 2D splines on \((r,\theta)\). 
* `max_iter` The maximum number of iterations possible for the batched CSR solver. 
* `res_tol` The residual tolerance for the batched CSR solver. Be careful! the relative residual provided here, will be used as "implicit residual" in ginkgo solver. 
* `batch_solver_logger` Indicates whether log information such as the residual and the number of iterations should be monitored. 
* `preconditioner_max_block_size` The maximum size of the Jacobi preconditioner used by the batched CSR solver.



**Template parameters:**


* `Mapping` A class describing a mapping from curvilinear coordinates to Cartesian coordinates. 




        

<hr>



### function operator() 

_Solve the Poisson-like equation._ 
```C++
template<class RHSFunction>
inline void PolarSplineFEMPoissonLikeSolver::operator() (
    RHSFunction const & rhs,
    PolarSplineRTheta spline
) const
```



This operator returns the coefficients associated with the B-Splines of the solution \(\phi\).




**Parameters:**


* `rhs` The rhs \(\rho\) of the Poisson-like equation. The type is templated but we can use the [**PoissonLikeRHSFunction**](classPoissonLikeRHSFunction.md) class. It must be an object with an operator() which evaluates a CoordRTheta and can be called from GPU. 
* `spline` The spline representation of the solution \(\phi\). 




        

<hr>



### function operator() 

_Solve the Poisson-like equation._ 
```C++
template<class RHSFunction>
inline void PolarSplineFEMPoissonLikeSolver::operator() (
    RHSFunction const & rhs,
    DFieldRTheta phi
) const
```



This operator uses the other operator () and returns the values on the grid of the solution \(\phi\).




**Parameters:**


* `rhs` The rhs \(\rho\) of the Poisson-like equation. The type is templated but we can use the [**PoissonLikeRHSFunction**](classPoissonLikeRHSFunction.md) class. It must be an object with an operator() which evaluates a CoordRTheta and can be called from GPU. 
* `phi` The values of the solution \(\phi\) on the given coords\_eval. 




        

<hr>



### function update\_coefficients 

_Update the coefficients_ \(alpha\) _and_\(beta\) _that define the equation._
```C++
inline void PolarSplineFEMPoissonLikeSolver::update_coefficients (
    ConstSpline2D coeff_alpha,
    ConstSpline2D coeff_beta
) 
```





**Parameters:**


* `coeff_alpha` The spline representation of the \(\alpha\) function in the definition of the Poisson-like equation. 
* `coeff_beta` The spline representation of the \(\beta\) function in the definition of the Poisson-like equation. 




        

<hr>
## Public Static Functions Documentation




### function get\_polar\_bspline\_vals 

_Get the value of the specified polar bspline at the specified point._ 
```C++
static inline KOKKOS_INLINE_FUNCTION double PolarSplineFEMPoissonLikeSolver::get_polar_bspline_vals (
    CoordRTheta coord,
    IdxBSPolar idx
) 
```





**Parameters:**


* `coord` The coordinate where the value of the polar bspline should be calculated. 
* `idx` The polar bspline of interest. 



**Returns:**

The value of the polar bspline at the coordinate. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/pde_solvers/polarpoissonlikesolver.hpp`

