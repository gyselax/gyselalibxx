

# Class PolarSplineFEMPoissonLikeSolver

**template &lt;class [**GridR**](structGridR.md), class [**GridTheta**](structGridTheta.md), class [**PolarBSplinesRTheta**](structPolarBSplinesRTheta.md), class SplineRThetaEvaluatorNullBound, class IdxRangeFull&gt;**



[**ClassList**](annotated.md) **>** [**PolarSplineFEMPoissonLikeSolver**](classPolarSplineFEMPoissonLikeSolver.md)



_Define a polar PDE solver for a Poisson-like equation._ [More...](#detailed-description)

* `#include <polarpoissonlikesolver.hpp>`















## Classes

| Type | Name |
| ---: | :--- |
| struct | [**EvalDeriv1DType**](structPolarSplineFEMPoissonLikeSolver_1_1EvalDeriv1DType.md) <br>_Object storing a value and a value of the derivative of a 1D function._  |
| struct | [**EvalDeriv2DType**](structPolarSplineFEMPoissonLikeSolver_1_1EvalDeriv2DType.md) <br>_Object storing a value and a value of the derivatives in each direction of a 2D function._  |
| struct | [**QDimRMesh**](structPolarSplineFEMPoissonLikeSolver_1_1QDimRMesh.md) <br>_Tag the first dimension for the quadrature mesh._  |
| struct | [**QDimThetaMesh**](structPolarSplineFEMPoissonLikeSolver_1_1QDimThetaMesh.md) <br>_Tag the second dimension for the quadrature mesh._  |
| struct | [**RBasisSubset**](structPolarSplineFEMPoissonLikeSolver_1_1RBasisSubset.md) <br> |
| struct | [**RCellDim**](structPolarSplineFEMPoissonLikeSolver_1_1RCellDim.md) <br> |
| struct | [**ThetaBasisSubset**](structPolarSplineFEMPoissonLikeSolver_1_1ThetaBasisSubset.md) <br> |
| struct | [**ThetaCellDim**](structPolarSplineFEMPoissonLikeSolver_1_1ThetaCellDim.md) <br> |


## Public Types

| Type | Name |
| ---: | :--- |
| typedef Idx&lt; [**RCellDim**](structPolarSplineFEMPoissonLikeSolver_1_1RCellDim.md), [**ThetaCellDim**](structPolarSplineFEMPoissonLikeSolver_1_1ThetaCellDim.md) &gt; | [**IdxCell**](#typedef-idxcell)  <br>_Tag an index of cell._  |
| typedef typename GridR::continuous\_dimension\_type | [**R**](#typedef-r)  <br>_The radial dimension._  |
| typedef typename GridTheta::continuous\_dimension\_type | [**Theta**](#typedef-theta)  <br>_The poloidal dimension._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**PolarSplineFEMPoissonLikeSolver**](#function-polarsplinefempoissonlikesolver) (ConstSpline2D coeff\_alpha, ConstSpline2D coeff\_beta, Mapping const & mapping, SplineRThetaEvaluatorNullBound const & spline\_evaluator, std::optional&lt; int &gt; max\_iter=std::nullopt, std::optional&lt; double &gt; res\_tol=std::nullopt, std::optional&lt; bool &gt; batch\_solver\_logger=std::nullopt, std::optional&lt; int &gt; preconditioner\_max\_block\_size=std::nullopt) <br>_Instantiate a polar Poisson-like solver using FEM with B-splines._  |
|  void | [**compute\_overlapping\_singular\_elements**](#function-compute_overlapping_singular_elements) (ConstSpline2D coeff\_alpha, ConstSpline2D coeff\_beta, Mapping const & mapping, SplineRThetaEvaluatorNullBound const & spline\_evaluator, Kokkos::View&lt; double \*\*, Kokkos::LayoutRight, Kokkos::HostSpace &gt; const values\_csr\_host, Kokkos::View&lt; int \*, Kokkos::LayoutRight, Kokkos::HostSpace &gt; const col\_idx\_csr\_host, Kokkos::View&lt; int \*, Kokkos::LayoutRight, Kokkos::HostSpace &gt; const nnz\_per\_row\_csr\_host) <br>_Computes the matrix element corresponding to singular elements overlapping with regular grid._  |
|  void | [**compute\_singular\_elements**](#function-compute_singular_elements) (ConstSpline2D coeff\_alpha, ConstSpline2D coeff\_beta, Mapping const & mapping, SplineRThetaEvaluatorNullBound const & spline\_evaluator, Kokkos::View&lt; double \*\*, Kokkos::LayoutRight, Kokkos::HostSpace &gt; const values\_csr\_host, Kokkos::View&lt; int \*, Kokkos::LayoutRight, Kokkos::HostSpace &gt; const col\_idx\_csr\_host, Kokkos::View&lt; int \*, Kokkos::LayoutRight, Kokkos::HostSpace &gt; const nnz\_per\_row\_csr\_host) <br>_Computes the matrix element corresponding to the singular area. ie: the region enclosing the O-point._  |
|  void | [**compute\_stencil\_elements**](#function-compute_stencil_elements) (ConstSpline2D coeff\_alpha, ConstSpline2D coeff\_beta, Mapping const & mapping, SplineRThetaEvaluatorNullBound const & spline\_evaluator, Kokkos::View&lt; double \*\*, Kokkos::LayoutRight, Kokkos::HostSpace &gt; const values\_csr\_host, Kokkos::View&lt; int \*, Kokkos::LayoutRight, Kokkos::HostSpace &gt; const col\_idx\_csr\_host, Kokkos::View&lt; int \*, Kokkos::LayoutRight, Kokkos::HostSpace &gt; const nnz\_per\_row\_csr\_host) <br>_Computes the matrix element corresponding to the regular stencil ie: out to singular or overlapping areas._  |
|  double | [**get\_matrix\_stencil\_element**](#function-get_matrix_stencil_element) (IdxBSRTheta idx\_test, IdxBSRTheta idx\_trial, ConstSpline2D coeff\_alpha, ConstSpline2D coeff\_beta, SplineRThetaEvaluatorNullBound const & evaluator, Mapping const & mapping) <br>_Computes the matrix element corresponding to two tensor product splines with index idx\_test and idx\_trial._  |
|  void | [**init\_nnz\_per\_line**](#function-init_nnz_per_line) (Kokkos::View&lt; int \*, Kokkos::LayoutRight &gt; nnz\_per\_row) const<br>_Fills the nnz data structure by computing the number of non-zero per line. This number is linked to the weak formulation and depends on_ \((r,\theta)\) _splines. After this function the array will contain: nnz\_per\_row[0] = 0. nnz\_per\_row[1] = 0. nnz\_per\_row[2] = number of non-zero elements in line 0. nnz\_per\_row[3] = number of non-zero elements in lines 0-1. ...\_per\_row nnz\_per\_row[matrix\_size] = number of non-zero elements in lines 0-(matrix\_size-1)._ |
|  void | [**operator()**](#function-operator) (RHSFunction const & rhs, host\_t&lt; PolarSplineRTheta &gt; spline) const<br>_Solve the Poisson-like equation._  |
|  void | [**operator()**](#function-operator_1) (RHSFunction const & rhs, DFieldRTheta phi) const<br>_Solve the Poisson-like equation._  |


## Public Static Functions

| Type | Name |
| ---: | :--- |
|  KOKKOS\_FUNCTION IdxRangeQuadratureRTheta | [**get\_quadrature\_points\_in\_cell**](#function-get_quadrature_points_in_cell) (int cell\_idx\_r, int cell\_idx\_theta) <br>_compute the quadrature range for a given pair of indices_  |
|  KOKKOS\_INLINE\_FUNCTION void | [**get\_value\_and\_gradient**](#function-get_value_and_gradient-12) (double & value, [**DVector**](classTensor.md)&lt; [**R\_cov**](structR__cov.md), [**Theta\_cov**](structTheta__cov.md) &gt; & derivs, [**EvalDeriv1DType**](structPolarSplineFEMPoissonLikeSolver_1_1EvalDeriv1DType.md) const & r\_basis, [**EvalDeriv1DType**](structPolarSplineFEMPoissonLikeSolver_1_1EvalDeriv1DType.md) const & theta\_basis) <br>_Computes the value and gradient from r\_basis and theta\_basis inputs._  |
|  KOKKOS\_INLINE\_FUNCTION void | [**get\_value\_and\_gradient**](#function-get_value_and_gradient-22) (double & value, [**DVector**](classTensor.md)&lt; [**R\_cov**](structR__cov.md), [**Theta\_cov**](structTheta__cov.md) &gt; & derivs, [**EvalDeriv2DType**](structPolarSplineFEMPoissonLikeSolver_1_1EvalDeriv2DType.md) const & basis, [**EvalDeriv2DType**](structPolarSplineFEMPoissonLikeSolver_1_1EvalDeriv2DType.md) const &) <br>_Computes the value and gradient from r\_basis and theta\_basis inputs._  |
|  KOKKOS\_FUNCTION double | [**templated\_weak\_integral\_element**](#function-templated_weak_integral_element) (IdxQuadratureR idx\_r, IdxQuadratureTheta idx\_theta, TestValDerivType const & test\_bspline\_val\_and\_deriv, TrialValDerivType const & trial\_bspline\_val\_and\_deriv, TestValDerivType const & test\_bspline\_val\_and\_deriv\_theta, TrialValDerivType const & trial\_bspline\_val\_and\_deriv\_theta, ConstSpline2D coeff\_alpha, ConstSpline2D coeff\_beta, SplineRThetaEvaluatorNullBound const & spline\_evaluator, Mapping const & mapping, DField&lt; IdxRangeQuadratureRTheta &gt; int\_volume) <br>_Computes a quadrature summand corresponding to the inner product._  |
|  KOKKOS\_FUNCTION int | [**theta\_mod**](#function-theta_mod) (int idx\_theta) <br>_Calculates the modulo idx\_theta in relation to cells number along_ \(\theta\) _direction ._ |
|  KOKKOS\_FUNCTION double | [**weak\_integral\_element**](#function-weak_integral_element) (IdxQuadratureR idx\_r, IdxQuadratureTheta idx\_theta, [**EvalDeriv2DType**](structPolarSplineFEMPoissonLikeSolver_1_1EvalDeriv2DType.md) const & test\_bspline\_val\_and\_deriv, [**EvalDeriv2DType**](structPolarSplineFEMPoissonLikeSolver_1_1EvalDeriv2DType.md) const & trial\_bspline\_val\_and\_deriv, ConstSpline2D coeff\_alpha, ConstSpline2D coeff\_beta, SplineRThetaEvaluatorNullBound const & evaluator, Mapping const & mapping, DField&lt; IdxRangeQuadratureRTheta &gt; int\_volume) <br>_compute the weak integral value._  |


























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
* `IdxRangeFull` The full index range of \(\phi\) including any batch dimensions. 




    
## Public Types Documentation




### typedef IdxCell 

_Tag an index of cell._ 
```C++
using PolarSplineFEMPoissonLikeSolver< GridR, GridTheta, PolarBSplinesRTheta, SplineRThetaEvaluatorNullBound, IdxRangeFull >::IdxCell =  Idx<RCellDim, ThetaCellDim>;
```




<hr>



### typedef R 

_The radial dimension._ 
```C++
using PolarSplineFEMPoissonLikeSolver< GridR, GridTheta, PolarBSplinesRTheta, SplineRThetaEvaluatorNullBound, IdxRangeFull >::R =  typename GridR::continuous_dimension_type;
```




<hr>



### typedef Theta 

_The poloidal dimension._ 
```C++
using PolarSplineFEMPoissonLikeSolver< GridR, GridTheta, PolarBSplinesRTheta, SplineRThetaEvaluatorNullBound, IdxRangeFull >::Theta =  typename GridTheta::continuous_dimension_type;
```




<hr>
## Public Functions Documentation




### function PolarSplineFEMPoissonLikeSolver 

_Instantiate a polar Poisson-like solver using FEM with B-splines._ 
```C++
template<class Mapping>
inline PolarSplineFEMPoissonLikeSolver::PolarSplineFEMPoissonLikeSolver (
    ConstSpline2D coeff_alpha,
    ConstSpline2D coeff_beta,
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


* `coeff_alpha` The spline representation of the \(\alpha\) function in the definition of the Poisson-like equation. 
* `coeff_beta` The spline representation of the \(\beta\) function in the definition of the Poisson-like equation. 
* `mapping` The mapping from the logical domain to the physical domain where the equation is defined. 
* `spline_evaluator` An evaluator for evaluating 2D splines on \((r,\theta)\). 
* `max_iter` The maximum number of iterations possible for the batched CSR solver. 
* `res_tol` The residual tolerance for the batched CSR solver. Be careful! the relative residual provided here, will be used as "implicit residual" in ginkgo solver. 
* `batch_solver_logger` Indicates whether log information such as the residual and the number of iterations should be monitored. 
* `preconditioner_max_block_size` The maximum size of the Jacobi preconditioner used by the batched CSR solver.



**Template parameters:**


* `Mapping` A class describing a mapping from curvilinear coordinates to Cartesian coordinates. 




        

<hr>



### function compute\_overlapping\_singular\_elements 

_Computes the matrix element corresponding to singular elements overlapping with regular grid._ 
```C++
template<class Mapping>
inline void PolarSplineFEMPoissonLikeSolver::compute_overlapping_singular_elements (
    ConstSpline2D coeff_alpha,
    ConstSpline2D coeff_beta,
    Mapping const & mapping,
    SplineRThetaEvaluatorNullBound const & spline_evaluator,
    Kokkos::View< double **, Kokkos::LayoutRight, Kokkos::HostSpace > const values_csr_host,
    Kokkos::View< int *, Kokkos::LayoutRight, Kokkos::HostSpace > const col_idx_csr_host,
    Kokkos::View< int *, Kokkos::LayoutRight, Kokkos::HostSpace > const nnz_per_row_csr_host
) 
```





**Parameters:**


* `coeff_alpha` The spline representation of the \(\alpha\) function in the definition of the Poisson-like equation. 
* `coeff_beta` The spline representation of the \(\beta\) function in the definition of the Poisson-like equation. 
* `mapping` The mapping from the logical domain to the physical domain where the equation is defined. 
* `spline_evaluator` An evaluator for evaluating 2D splines on \((r,\theta)\). 
* `values_csr_host` A 2D Kokkos view which stores the values of non-zero elements for the whole batch. 
* `col_idx_csr_host` A 1D Kokkos view which stores the column indices for each non-zero component.(only for one matrix) 
* `nnz_per_row_csr_host` A 1D Kokkos view of length matrix\_size+1 which stores the count of the non-zeros along the lines of the matrix. 




        

<hr>



### function compute\_singular\_elements 

_Computes the matrix element corresponding to the singular area. ie: the region enclosing the O-point._ 
```C++
template<class Mapping>
inline void PolarSplineFEMPoissonLikeSolver::compute_singular_elements (
    ConstSpline2D coeff_alpha,
    ConstSpline2D coeff_beta,
    Mapping const & mapping,
    SplineRThetaEvaluatorNullBound const & spline_evaluator,
    Kokkos::View< double **, Kokkos::LayoutRight, Kokkos::HostSpace > const values_csr_host,
    Kokkos::View< int *, Kokkos::LayoutRight, Kokkos::HostSpace > const col_idx_csr_host,
    Kokkos::View< int *, Kokkos::LayoutRight, Kokkos::HostSpace > const nnz_per_row_csr_host
) 
```





**Parameters:**


* `coeff_alpha` The spline representation of the \(\alpha\) function in the definition of the Poisson-like equation. 
* `coeff_beta` The spline representation of the \(\beta\) function in the definition of the Poisson-like equation. 
* `mapping` The mapping from the logical domain to the physical domain where the equation is defined. 
* `spline_evaluator` An evaluator for evaluating 2D splines on \((r,\theta)\). 
* `values_csr_host` A 2D Kokkos view which stores the values of non-zero elements for the whole batch. 
* `col_idx_csr_host` A 1D Kokkos view which stores the column indices for each non-zero component.(only for one matrix). 
* `nnz_per_row_csr_host` A 1D Kokkos view of length matrix\_size+1 which stores the count of the non-zeros along the lines of the matrix. 




        

<hr>



### function compute\_stencil\_elements 

_Computes the matrix element corresponding to the regular stencil ie: out to singular or overlapping areas._ 
```C++
template<class Mapping>
inline void PolarSplineFEMPoissonLikeSolver::compute_stencil_elements (
    ConstSpline2D coeff_alpha,
    ConstSpline2D coeff_beta,
    Mapping const & mapping,
    SplineRThetaEvaluatorNullBound const & spline_evaluator,
    Kokkos::View< double **, Kokkos::LayoutRight, Kokkos::HostSpace > const values_csr_host,
    Kokkos::View< int *, Kokkos::LayoutRight, Kokkos::HostSpace > const col_idx_csr_host,
    Kokkos::View< int *, Kokkos::LayoutRight, Kokkos::HostSpace > const nnz_per_row_csr_host
) 
```





**Parameters:**


* `coeff_alpha` The spline representation of the \(\alpha\) function in the definition of the Poisson-like equation. 
* `coeff_beta` The spline representation of the \(\beta\) function in the definition of the Poisson-like equation. 
* `mapping` The mapping from the logical domain to the physical domain where the equation is defined. 
* `spline_evaluator` An evaluator for evaluating 2D splines on \((r,\theta)\). 
* `values_csr_host` A 2D Kokkos view which stores the values of non-zero elements for the whole batch. 
* `col_idx_csr_host` A 1D Kokkos view which stores the column indices for each non-zero component.(only for one matrix) 
* `nnz_per_row_csr_host` A 1D Kokkos view of length matrix\_size+1 which stores the count of the non-zeros along the lines of the matrix. 




        

<hr>



### function get\_matrix\_stencil\_element 

_Computes the matrix element corresponding to two tensor product splines with index idx\_test and idx\_trial._ 
```C++
template<class Mapping>
inline double PolarSplineFEMPoissonLikeSolver::get_matrix_stencil_element (
    IdxBSRTheta idx_test,
    IdxBSRTheta idx_trial,
    ConstSpline2D coeff_alpha,
    ConstSpline2D coeff_beta,
    SplineRThetaEvaluatorNullBound const & evaluator,
    Mapping const & mapping
) 
```





**Parameters:**


* `idx_test` The index for polar B-spline in the test space. 
* `idx_trial` The index for polar B-spline in the trial space. 
* `coeff_alpha` The spline representation of the \(\alpha\) function in the definition of the Poisson-like equation. 
* `coeff_beta` The spline representation of the \(\beta\) function in the definition of the Poisson-like equation. 
* `evaluator` An evaluator for evaluating 2D splines on \((r, \theta)\). 
* `mapping` The mapping from the logical domain to the physical domain where the equation is defined. 



**Returns:**

The value of the matrix element. 





        

<hr>



### function init\_nnz\_per\_line 

_Fills the nnz data structure by computing the number of non-zero per line. This number is linked to the weak formulation and depends on_ \((r,\theta)\) _splines. After this function the array will contain: nnz\_per\_row[0] = 0. nnz\_per\_row[1] = 0. nnz\_per\_row[2] = number of non-zero elements in line 0. nnz\_per\_row[3] = number of non-zero elements in lines 0-1. ...\_per\_row nnz\_per\_row[matrix\_size] = number of non-zero elements in lines 0-(matrix\_size-1)._
```C++
inline void PolarSplineFEMPoissonLikeSolver::init_nnz_per_line (
    Kokkos::View< int *, Kokkos::LayoutRight > nnz_per_row
) const
```





**Parameters:**


* `nnz_per_row` A 1D Kokkos view of length matrix\_size+1 which stores the sum of the non-zeros in the matrix on all lines up to the one in. 




        

<hr>



### function operator() 

_Solve the Poisson-like equation._ 
```C++
template<class RHSFunction>
inline void PolarSplineFEMPoissonLikeSolver::operator() (
    RHSFunction const & rhs,
    host_t< PolarSplineRTheta > spline
) const
```



This operator returns the coefficients associated with the B-Splines of the solution \(\phi\).




**Parameters:**


* `rhs` The rhs \(\rho\) of the Poisson-like equation. The type is templated but we can use the [**PoissonLikeRHSFunction**](classPoissonLikeRHSFunction.md) class. It must be an object with an operator() which evaluates a CoordRTheta and can be called from CPU. 
* `spline` The spline representation of the solution \(\phi\), also used as initial data for the iterative solver. 




        

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


* `rhs` The rhs \(\rho\) of the Poisson-like equation. The type is templated but we can use the [**PoissonLikeRHSFunction**](classPoissonLikeRHSFunction.md) class. It must be an object with an operator() which evaluates a CoordRTheta and can be called from CPU. 
* `phi` The values of the solution \(\phi\) on the given coords\_eval, also used as initial data for the iterative solver. 




        

<hr>
## Public Static Functions Documentation




### function get\_quadrature\_points\_in\_cell 

_compute the quadrature range for a given pair of indices_ 
```C++
static inline KOKKOS_FUNCTION IdxRangeQuadratureRTheta PolarSplineFEMPoissonLikeSolver::get_quadrature_points_in_cell (
    int cell_idx_r,
    int cell_idx_theta
) 
```





**Parameters:**


* `cell_idx_r` The index for radial direction 
* `cell_idx_theta` The index for poloidal direction 



**Returns:**

The quadrature range corresponding to the \((r,\theta)\) indices. 





        

<hr>



### function get\_value\_and\_gradient [1/2]

_Computes the value and gradient from r\_basis and theta\_basis inputs._ 
```C++
static inline KOKKOS_INLINE_FUNCTION void PolarSplineFEMPoissonLikeSolver::get_value_and_gradient (
    double & value,
    DVector < R_cov , Theta_cov > & derivs,
    EvalDeriv1DType const & r_basis,
    EvalDeriv1DType const & theta_basis
) 
```





**Parameters:**


* `value` The product of radial and poloidal values.
* `derivs` derivatives over \((r, \theta)\) directions.
* `r_basis` A data structure containing values and derivative over radial direction.
* `theta_basis` A data structure containing values and derivative over poloidal direction. 




        

<hr>



### function get\_value\_and\_gradient [2/2]

_Computes the value and gradient from r\_basis and theta\_basis inputs._ 
```C++
static inline KOKKOS_INLINE_FUNCTION void PolarSplineFEMPoissonLikeSolver::get_value_and_gradient (
    double & value,
    DVector < R_cov , Theta_cov > & derivs,
    EvalDeriv2DType const & basis,
    EvalDeriv2DType const &
) 
```





**Parameters:**


* `value` The product of radial and poloidal values.
* `derivs` derivatives over \((r, \theta)\) directions.
* `basis` A data structure containing values and derivative over radial and poloidal directions. 




        

<hr>



### function templated\_weak\_integral\_element 

_Computes a quadrature summand corresponding to the inner product._ 
```C++
template<class Mapping, class TestValDerivType, class TrialValDerivType>
static inline KOKKOS_FUNCTION double PolarSplineFEMPoissonLikeSolver::templated_weak_integral_element (
    IdxQuadratureR idx_r,
    IdxQuadratureTheta idx_theta,
    TestValDerivType const & test_bspline_val_and_deriv,
    TrialValDerivType const & trial_bspline_val_and_deriv,
    TestValDerivType const & test_bspline_val_and_deriv_theta,
    TrialValDerivType const & trial_bspline_val_and_deriv_theta,
    ConstSpline2D coeff_alpha,
    ConstSpline2D coeff_beta,
    SplineRThetaEvaluatorNullBound const & spline_evaluator,
    Mapping const & mapping,
    DField< IdxRangeQuadratureRTheta > int_volume
) 
```





**Parameters:**


* `idx_r` The index for radial direction. 
* `idx_theta` The index for poloidal direction 
* `test_bspline_val_and_deriv` The data structure containing the derivatives over radial and poloidal directions for test space. 
* `trial_bspline_val_and_deriv` The data structure containing the derivatives over radial and poloidal directions for trial space. 
* `test_bspline_val_and_deriv_theta` The data structure containing the value and derivative along poloidal direction for test space. 
* `trial_bspline_val_and_deriv_theta` The data structure containing the value and derivative along poloidal direction for trial space. 
* `coeff_alpha` The spline representation of the \(\alpha\) function in the definition of the Poisson-like equation. 
* `coeff_beta` The spline representation of the \(\beta\) function in the definition of the Poisson-like equation. 
* `spline_evaluator` An evaluator for evaluating 2D splines on \((r,\theta)\). 
* `mapping` The mapping from the logical domain to the physical domain where the equation is defined. 
* `int_volume` The integral volume associated with each point used in the quadrature scheme. 



**Returns:**

inner product of the test and trial spline is computed using a quadrature. This function returns one summand of the quadrature for the quadrature point given by the indices. 





        

<hr>



### function theta\_mod 

_Calculates the modulo idx\_theta in relation to cells number along_ \(\theta\) _direction ._
```C++
static inline KOKKOS_FUNCTION int PolarSplineFEMPoissonLikeSolver::theta_mod (
    int idx_theta
) 
```





**Parameters:**


* `idx_theta` \(\theta\) index.



**Returns:**

The corresponding indice modulo \(\theta\) direction cells number 





        

<hr>



### function weak\_integral\_element 

_compute the weak integral value._ 
```C++
template<class Mapping>
static inline KOKKOS_FUNCTION double PolarSplineFEMPoissonLikeSolver::weak_integral_element (
    IdxQuadratureR idx_r,
    IdxQuadratureTheta idx_theta,
    EvalDeriv2DType const & test_bspline_val_and_deriv,
    EvalDeriv2DType const & trial_bspline_val_and_deriv,
    ConstSpline2D coeff_alpha,
    ConstSpline2D coeff_beta,
    SplineRThetaEvaluatorNullBound const & evaluator,
    Mapping const & mapping,
    DField< IdxRangeQuadratureRTheta > int_volume
) 
```





**Parameters:**


* `idx_r` The index for radial direction. 
* `idx_theta` The index for poloidal direction 
* `test_bspline_val_and_deriv` The data structure containing the derivatives over radial and poloidal directions for test space. 
* `trial_bspline_val_and_deriv` The data structure containing the derivatives over radial and poloidal directions for trial space. 
* `coeff_alpha` The spline representation of the \(\alpha\) function in the definition of the Poisson-like equation. 
* `coeff_beta` The spline representation of the \(\beta\) function in the definition of the Poisson-like equation. 
* `mapping` The mapping from the logical domain to the physical domain where the equation is defined. 
* `evaluator` An evaluator for evaluating 2D splines on \((r,\theta)\). 
* `int_volume` The integral volume associated with each point used in the quadrature scheme. 



**Returns:**

The value of the weak integral. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/pde_solvers/polarpoissonlikesolver.hpp`

