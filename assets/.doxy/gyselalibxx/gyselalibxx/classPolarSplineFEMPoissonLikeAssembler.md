

# Class PolarSplineFEMPoissonLikeAssembler

**template &lt;typename [**GridR**](structGridR.md), typename [**GridTheta**](structGridTheta.md), typename [**PolarBSplinesRTheta**](structPolarBSplinesRTheta.md), typename SplineRThetaEvaluatorNullBound, typename QDimRMesh, typename QDimThetaMesh, class IdxRangeFull&gt;**



[**ClassList**](annotated.md) **>** [**PolarSplineFEMPoissonLikeAssembler**](classPolarSplineFEMPoissonLikeAssembler.md)



_An operator to assemble a Poisson-like stiffness matrix using polar B-splines._ [More...](#detailed-description)

* `#include <polarpoissonlikeassembler.hpp>`















## Classes

| Type | Name |
| ---: | :--- |
| struct | [**InternalBatchDim**](structPolarSplineFEMPoissonLikeAssembler_1_1InternalBatchDim.md) <br>_The tag for the batch dimension for the equation. This is public due to Cuda._  |


## Public Types

| Type | Name |
| ---: | :--- |
| typedef typename GridR::continuous\_dimension\_type | [**R**](#typedef-r)  <br>_The radial dimension._  |
| typedef typename GridTheta::continuous\_dimension\_type | [**Theta**](#typedef-theta)  <br>_The poloidal dimension._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**PolarSplineFEMPoissonLikeAssembler**](#function-polarsplinefempoissonlikeassembler) (Field&lt; double, IdxRangeQuadratureRTheta &gt; int\_volume) <br>_Instantiate the assembler operator._  |
|  void | [**compute\_overlapping\_singular\_elements**](#function-compute_overlapping_singular_elements) (ConstSpline2D coeff\_alpha, ConstSpline2D coeff\_beta, Mapping const & mapping, SplineRThetaEvaluatorNullBound const & spline\_evaluator, Kokkos::View&lt; double \*\*, Kokkos::LayoutRight, Kokkos::HostSpace &gt; const values\_csr\_host, Kokkos::View&lt; int \*, Kokkos::LayoutRight, Kokkos::HostSpace &gt; const col\_idx\_csr\_host, Kokkos::View&lt; int \*, Kokkos::LayoutRight, Kokkos::HostSpace &gt; const nnz\_per\_row\_csr\_host) <br>_Computes the matrix element corresponding to singular elements overlapping with regular grid._  |
|  void | [**compute\_singular\_elements**](#function-compute_singular_elements) (ConstSpline2D coeff\_alpha, ConstSpline2D coeff\_beta, Mapping const & mapping, SplineRThetaEvaluatorNullBound const & spline\_evaluator, Kokkos::View&lt; double \*\*, Kokkos::LayoutRight, Kokkos::HostSpace &gt; const values\_csr\_host, Kokkos::View&lt; int \*, Kokkos::LayoutRight, Kokkos::HostSpace &gt; const col\_idx\_csr\_host, Kokkos::View&lt; int \*, Kokkos::LayoutRight, Kokkos::HostSpace &gt; const nnz\_per\_row\_csr\_host) <br>_Computes the matrix element corresponding to the singular area. ie: the region enclosing the O-point._  |
|  void | [**compute\_stencil\_elements**](#function-compute_stencil_elements) (ConstSpline2D coeff\_alpha, ConstSpline2D coeff\_beta, Mapping const & mapping, SplineRThetaEvaluatorNullBound const & spline\_evaluator, Kokkos::View&lt; double \*\*, Kokkos::LayoutRight, Kokkos::HostSpace &gt; const values\_csr\_host, Kokkos::View&lt; int \*, Kokkos::LayoutRight, Kokkos::HostSpace &gt; const col\_idx\_csr\_host, Kokkos::View&lt; int \*, Kokkos::LayoutRight, Kokkos::HostSpace &gt; const nnz\_per\_row\_csr\_host) <br>_Computes the matrix element corresponding to the regular stencil ie: out to singular or overlapping areas._  |
|  double | [**get\_matrix\_stencil\_element**](#function-get_matrix_stencil_element) (IdxBSRTheta idx\_test, IdxBSRTheta idx\_trial, ConstSpline2D coeff\_alpha, ConstSpline2D coeff\_beta, SplineRThetaEvaluatorNullBound const & evaluator, Mapping const & mapping) <br>_Computes the matrix element corresponding to two tensor product splines with index idx\_test and idx\_trial._  |
|  void | [**init\_nnz\_per\_line**](#function-init_nnz_per_line) (Kokkos::View&lt; int \*, Kokkos::LayoutRight &gt; nnz\_per\_row) const<br>_Fills the nnz data structure by computing the number of non-zero per line. This number is linked to the weak formulation and depends on_ \((r,\theta)\) _splines. After this function the array will contain: nnz\_per\_row[0] = 0. nnz\_per\_row[1] = 0. nnz\_per\_row[2] = number of non-zero elements in line 0. nnz\_per\_row[3] = number of non-zero elements in lines 0-1. ...\_per\_row nnz\_per\_row[matrix\_size] = number of non-zero elements in lines 0-(matrix\_size-1)._ |
|  void | [**operator()**](#function-operator) (std::unique\_ptr&lt; [**MatrixBatchCsr**](classMatrixBatchCsr.md)&lt; Kokkos::DefaultExecutionSpace, MatrixBatchCsrSolver::CG &gt; &gt; & gko\_matrix, ConstSpline2D coeff\_alpha, ConstSpline2D coeff\_beta, Mapping const & mapping, SplineRThetaEvaluatorNullBound const & spline\_evaluator, std::optional&lt; int &gt; max\_iter=std::nullopt, std::optional&lt; double &gt; res\_tol=std::nullopt, std::optional&lt; bool &gt; batch\_solver\_logger=std::nullopt, std::optional&lt; int &gt; preconditioner\_max\_block\_size=std::nullopt) <br>_Assemble the stiffness matrix._  |


## Public Static Functions

| Type | Name |
| ---: | :--- |
|  KOKKOS\_INLINE\_FUNCTION IdxBSPolar | [**to\_polar**](#function-to_polar) (IdxBSRTheta idx) <br>_Convert a 2D (r,theta) bspline index into a polar bspline index._  |
|  KOKKOS\_FUNCTION double | [**weak\_integral\_element**](#function-weak_integral_element) (IdxBSPolar idx\_test, IdxBSPolar idx\_trial, IdxQuadratureRTheta idx\_quad, ConstSpline2D coeff\_alpha, ConstSpline2D coeff\_beta, SplineRThetaEvaluatorNullBound const & spline\_evaluator, Mapping const & mapping, DField&lt; IdxRangeQuadratureRTheta &gt; int\_volume) <br>_Computes a quadrature summand for computing the integral corresponding to the inner product._  |


























## Detailed Description


Assemble the finite element stiffness matrix with entries \(\int \alpha \nabla B_i \cdot \nabla B_j + \beta B_i B_j dx\) for \(0 < i,j < N\) where \(N\) is the number of polar splines. We only consider splines that vanish on the boundary.


For more details see the `PolarSplineFEMPoissonLikeSolver` and Emily Bourne's thesis ("Non-Uniform Numerical Schemes for the Modelling of Turbulence
in the 5D GYSELA Code". December 2022.) 

**Template parameters:**


* [**GridR**](structGridR.md) The radial grid type. 
* [**GridR**](structGridR.md) The poloidal grid type. 
* [**PolarBSplinesRTheta**](structPolarBSplinesRTheta.md) The type of the 2D polar B-splines (on the coordinate system \((r,\theta)\) including B-splines which traverse the O point). 
* `SplineRThetaEvaluatorNullBound` The type of the 2D (cross-product) spline evaluator. 
* `QDimRMesh` The radial quadrature grid type. 
* `QDimThetaMesh` The poloidal quadrature grid type. 
* `IdxRangeFull` The full index range of \(\phi\) including any batch dimensions. 




    
## Public Types Documentation




### typedef R 

_The radial dimension._ 
```C++
using PolarSplineFEMPoissonLikeAssembler< GridR, GridTheta, PolarBSplinesRTheta, SplineRThetaEvaluatorNullBound, QDimRMesh, QDimThetaMesh, IdxRangeFull >::R =  typename GridR::continuous_dimension_type;
```




<hr>



### typedef Theta 

_The poloidal dimension._ 
```C++
using PolarSplineFEMPoissonLikeAssembler< GridR, GridTheta, PolarBSplinesRTheta, SplineRThetaEvaluatorNullBound, QDimRMesh, QDimThetaMesh, IdxRangeFull >::Theta =  typename GridTheta::continuous_dimension_type;
```




<hr>
## Public Functions Documentation




### function PolarSplineFEMPoissonLikeAssembler 

_Instantiate the assembler operator._ 
```C++
inline explicit PolarSplineFEMPoissonLikeAssembler::PolarSplineFEMPoissonLikeAssembler (
    Field< double, IdxRangeQuadratureRTheta > int_volume
) 
```





**Parameters:**


* `int_volume` The initialised field of Jacobian values of the mapping. 




        

<hr>



### function compute\_overlapping\_singular\_elements 

_Computes the matrix element corresponding to singular elements overlapping with regular grid._ 
```C++
template<class Mapping>
inline void PolarSplineFEMPoissonLikeAssembler::compute_overlapping_singular_elements (
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
inline void PolarSplineFEMPoissonLikeAssembler::compute_singular_elements (
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
inline void PolarSplineFEMPoissonLikeAssembler::compute_stencil_elements (
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
inline double PolarSplineFEMPoissonLikeAssembler::get_matrix_stencil_element (
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
inline void PolarSplineFEMPoissonLikeAssembler::init_nnz_per_line (
    Kokkos::View< int *, Kokkos::LayoutRight > nnz_per_row
) const
```





**Parameters:**


* `nnz_per_row` A 1D Kokkos view of length matrix\_size+1 which stores the sum of the non-zeros in the matrix on all lines up to the one in. 




        

<hr>



### function operator() 

_Assemble the stiffness matrix._ 
```C++
template<typename Mapping>
inline void PolarSplineFEMPoissonLikeAssembler::operator() (
    std::unique_ptr< MatrixBatchCsr < Kokkos::DefaultExecutionSpace, MatrixBatchCsrSolver::CG > > & gko_matrix,
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





**Parameters:**


* `gko_matrix` The pointer to the assembled matrix. 
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
## Public Static Functions Documentation




### function to\_polar 

_Convert a 2D (r,theta) bspline index into a polar bspline index._ 
```C++
static inline KOKKOS_INLINE_FUNCTION IdxBSPolar PolarSplineFEMPoissonLikeAssembler::to_polar (
    IdxBSRTheta idx
) 
```





**Parameters:**


* `idx` The 2D (r,theta) bspline index. 



**Returns:**

The polar bspline index. 





        

<hr>



### function weak\_integral\_element 

_Computes a quadrature summand for computing the integral corresponding to the inner product._ 
```C++
template<class Mapping>
static inline KOKKOS_FUNCTION double PolarSplineFEMPoissonLikeAssembler::weak_integral_element (
    IdxBSPolar idx_test,
    IdxBSPolar idx_trial,
    IdxQuadratureRTheta idx_quad,
    ConstSpline2D coeff_alpha,
    ConstSpline2D coeff_beta,
    SplineRThetaEvaluatorNullBound const & spline_evaluator,
    Mapping const & mapping,
    DField< IdxRangeQuadratureRTheta > int_volume
) 
```



Inner product of the test and trial spline is computed using a quadrature. This function returns one summand of the quadrature for the quadrature point given by the indices.




**Parameters:**


* `idx_test` The index of the test basis spline. 
* `idx_trial` The index of the trial basis spline. 
* `idx_quad` The index for the point in the quadrature scheme. 
* `coeff_alpha` The spline representation of the \(\alpha\) function in the definition of the Poisson-like equation. 
* `coeff_beta` The spline representation of the \(\beta\) function in the definition of the Poisson-like equation. 
* `spline_evaluator` An evaluator for evaluating 2D splines on \((r,\theta)\). 
* `mapping` The mapping from the logical domain to the physical domain where the equation is defined. 
* `int_volume` The integral volume associated with each point used in the quadrature scheme. 



**Returns:**

Value of the quadrature summand 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/pde_solvers/polarpoissonlikeassembler.hpp`

