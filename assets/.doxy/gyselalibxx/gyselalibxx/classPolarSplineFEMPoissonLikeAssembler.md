

# Class PolarSplineFEMPoissonLikeAssembler

**template &lt;typename [**GridR**](structGridR.md), typename [**GridTheta**](structGridTheta.md), typename [**PolarBSplinesRTheta**](structPolarBSplinesRTheta.md), typename QDimRMesh, typename QDimThetaMesh, concepts::MetricTensor MetricTensorType, typename IdxRangeFull&gt;**



[**ClassList**](annotated.md) **>** [**PolarSplineFEMPoissonLikeAssembler**](classPolarSplineFEMPoissonLikeAssembler.md)



_An operator to assemble a Poisson-like stiffness matrix using polar B-splines._ [More...](#detailed-description)

* `#include <polar_spline_fem_poisson_like_assembler.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef typename GridR::continuous\_dimension\_type | [**R**](#typedef-r)  <br>_The radial dimension._  |
| typedef detail\_poisson:: [**TensorBSplineStash**](classPolarSplineFEMPoissonLikeAssembler.md#typedef-tensorbsplinestash)&lt; QDimRMesh, QDimThetaMesh, [**BSplinesR**](structBSplinesR.md), [**BSplinesTheta**](structBSplinesTheta.md) &gt; | [**TensorBSplineStash**](#typedef-tensorbsplinestash)  <br>_The values and derivatives of the 1D B-splines at the quadrature points._  |
| typedef typename GridTheta::continuous\_dimension\_type | [**Theta**](#typedef-theta)  <br>_The poloidal dimension._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**PolarSplineFEMPoissonLikeAssembler**](#function-polarsplinefempoissonlikeassembler) (Field&lt; double, IdxRangeQuadratureRTheta &gt; int\_volume, MetricTensorType const & metric\_tensor, IdxRangeFull idx\_range\_full) <br>_Instantiate the assembler operator._  |
|  void | [**compute\_singular\_singular\_col\_idx**](#function-compute_singular_singular_col_idx) (Kokkos::View&lt; int \*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace &gt; const col\_idx\_csr, Kokkos::View&lt; int \*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace &gt; const nnz\_per\_row\_csr) <br>_Computes the column indices of the matrix elements corresponding to the inner products of singular basis functions and singular basis functions._  |
|  void | [**compute\_singular\_singular\_elements**](#function-compute_singular_singular_elements) (CoeffAlpha const & coeff\_alpha, CoeffBeta const & coeff\_beta, MetricTensorType const & metric\_tensor, Kokkos::View&lt; double \*\*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace &gt; const values\_csr, Kokkos::View&lt; int \*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace &gt; const col\_idx\_csr, Kokkos::View&lt; int \*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace &gt; const nnz\_per\_row\_csr) <br>_Computes the matrix elements corresponding to the inner products of singular basis functions and singular basis functions._  |
|  void | [**compute\_singular\_tensor\_col\_idx**](#function-compute_singular_tensor_col_idx) (Kokkos::View&lt; int \*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace &gt; const col\_idx\_csr, Kokkos::View&lt; int \*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace &gt; const nnz\_per\_row\_csr) <br>_Computes the column indices of the matrix elements corresponding to the inner products of singular basis functions and tensor basis functions._  |
|  void | [**compute\_singular\_tensor\_elements**](#function-compute_singular_tensor_elements) (CoeffAlpha const & coeff\_alpha, CoeffBeta const & coeff\_beta, MetricTensorType const & metric\_tensor, Kokkos::View&lt; double \*\*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace &gt; const values\_csr, Kokkos::View&lt; int \*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace &gt; const col\_idx\_csr, Kokkos::View&lt; int \*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace &gt; const nnz\_per\_row\_csr) <br>_Computes the matrix elements corresponding to the inner products of singular basis functions and tensor basis functions._  |
|  void | [**compute\_tensor\_tensor\_col\_idx**](#function-compute_tensor_tensor_col_idx) (Kokkos::View&lt; int \*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace &gt; const col\_idx\_csr, Kokkos::View&lt; int \*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace &gt; const nnz\_per\_row\_csr) <br>_Computes the column indices of the matrix elements corresponding to the inner products of tensor basis functions and tensor basis functions._  |
|  void | [**compute\_tensor\_tensor\_elements**](#function-compute_tensor_tensor_elements) (CoeffAlpha const & coeff\_alpha, CoeffBeta const & coeff\_beta, MetricTensorType const & metric\_tensor, Kokkos::View&lt; double \*\*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace &gt; const values\_csr, Kokkos::View&lt; int \*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace &gt; const col\_idx\_csr, Kokkos::View&lt; int \*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace &gt; const nnz\_per\_row\_csr) <br>_Computes the matrix elements corresponding to the inner products of tensor basis functions and tensor basis functions._  |
|  [**TensorBSplineStash**](classPolarSplineFEMPoissonLikeAssembler.md#typedef-tensorbsplinestash) | [**get\_tensor\_bspline\_stash**](#function-get_tensor_bspline_stash) () const<br>_Get the_ [_**TensorBSplineStash**_](classPolarSplineFEMPoissonLikeAssembler.md#typedef-tensorbsplinestash) _of the values and derivatives of the 1D B-splines at the quadrature points._ |
|  void | [**init\_nnz\_per\_line**](#function-init_nnz_per_line) (Kokkos::View&lt; int \*, Kokkos::LayoutRight &gt; nnz\_per\_row) const<br>_Fills the nnz data structure by computing the number of non-zero per line. This number is linked to the weak formulation and depends on_ \((r,\theta)\) _splines. After this function the array will contain: nnz\_per\_row[0] = 0. nnz\_per\_row[1] = 0. nnz\_per\_row[2] = number of non-zero elements in line 0. nnz\_per\_row[3] = number of non-zero elements in lines 0-1. ...\_per\_row nnz\_per\_row[matrix\_size] = number of non-zero elements in lines 0-(matrix\_size-1)._ |
|  void | [**operator()**](#function-operator) (std::unique\_ptr&lt; [**MatrixBatchCsr**](classMatrixBatchCsr.md)&lt; Kokkos::DefaultExecutionSpace, MatrixBatchCsrSolver::CG &gt; &gt; const & gko\_matrix, CoeffAlpha const & coeff\_alpha, CoeffBeta const & coeff\_beta) <br>_Assemble the stiffness matrix._  |
|  void | [**setup\_sparse\_matrix**](#function-setup_sparse_matrix) (std::unique\_ptr&lt; [**MatrixBatchCsr**](classMatrixBatchCsr.md)&lt; Kokkos::DefaultExecutionSpace, MatrixBatchCsrSolver::CG &gt; &gt; & gko\_matrix, std::optional&lt; int &gt; max\_iter=std::nullopt, std::optional&lt; double &gt; res\_tol=std::nullopt, std::optional&lt; bool &gt; batch\_solver\_logger=std::nullopt, std::optional&lt; int &gt; preconditioner\_max\_block\_size=std::nullopt, std::optional&lt; GkoMatrixResidualNorm &gt; residual\_norm=std::nullopt) <br>_Compute the sparsity pattern for the stiffness matrix._  |


## Public Static Functions

| Type | Name |
| ---: | :--- |
|  KOKKOS\_FUNCTION void | [**fill\_test\_weights**](#function-fill_test_weights) (const Kokkos::TeamPolicy&lt;&gt;::member\_type & team, QuadField test\_mass, DVectorField&lt; IdxRangeQuadratureRTheta, VectorIndexSet&lt; [**R**](classPolarSplineFEMPoissonLikeAssembler.md#typedef-r)[**R**](structR.md), [**Theta**](classPolarSplineFEMPoissonLikeAssembler.md#typedef-theta)[**Theta**](structTheta.md) &gt; &gt; test\_stiffness, IdxBSRTheta idx\_test, IdxBatch idx\_batch, CoeffAlpha const & coeff\_alpha, CoeffBeta const & coeff\_beta, MetricTensorType const & metric\_tensor, IdxRangeQuadratureRTheta const & full\_quad\_idx\_range, DField&lt; IdxRangeQuadratureRTheta &gt; int\_volume, [**TensorBSplineStash**](classPolarSplineFEMPoissonLikeAssembler.md#typedef-tensorbsplinestash) const & bspline\_stash) <br>_Computes the test\_mass\_weight and test\_stiffness\_weight at every quadrature point in the support of a test function and saves them in team scratch memory._  |
|  KOKKOS\_FUNCTION double | [**get\_matrix\_stencil\_element**](#function-get_matrix_stencil_element) (const Kokkos::TeamPolicy&lt;&gt;::member\_type & team, IdxBSRTheta idx\_test, IdxBSRTheta idx\_trial, QuadField test\_mass, DVectorField&lt; IdxRangeQuadratureRTheta, VectorIndexSet&lt; [**R**](classPolarSplineFEMPoissonLikeAssembler.md#typedef-r)[**R**](structR.md), [**Theta**](classPolarSplineFEMPoissonLikeAssembler.md#typedef-theta)[**Theta**](structTheta.md) &gt; &gt; test\_stiffness, IdxRangeQuadratureRTheta const & full\_quad\_idx\_range, [**TensorBSplineStash**](classPolarSplineFEMPoissonLikeAssembler.md#typedef-tensorbsplinestash) const & bspline\_stash) <br>_Computes the matrix element corresponding to two tensor product splines with index idx\_test and idx\_trial._  |
|  KOKKOS\_FUNCTION std::tuple&lt; double, DVector&lt; [**R**](classPolarSplineFEMPoissonLikeAssembler.md#typedef-r)[**R**](structR.md), [**Theta**](classPolarSplineFEMPoissonLikeAssembler.md#typedef-theta)[**Theta**](structTheta.md) &gt; &gt; | [**get\_test\_weights**](#function-get_test_weights) (IdxTestType idx\_test, IdxQuadratureRTheta idx\_quad, IdxBatch idx\_batch, CoeffAlpha const & coeff\_alpha, CoeffBeta const & coeff\_beta, MetricTensorType const & metric\_tensor, DField&lt; IdxRangeQuadratureRTheta &gt; int\_volume, [**TensorBSplineStash**](classPolarSplineFEMPoissonLikeAssembler.md#typedef-tensorbsplinestash) const & bspline\_stash) <br>_Computes the test\_stiffness\_weight and test\_mass\_weight for a given quadrature point._  |
|  void | [**stash\_basis\_eval**](#function-stash_basis_eval) (DField&lt; IdxRange&lt; QDimMesh, BSplines, DerivDim &gt; &gt; basis\_vals, Field&lt; Idx&lt; BSplines &gt;, IdxRange&lt; QDimMesh &gt; &gt; first\_basis) <br>_Fill the stash with the values and derivatives of the non-zero 1D B-splines at the 1D quadrature points._  |
|  KOKKOS\_INLINE\_FUNCTION IdxBSPolar | [**to\_polar**](#function-to_polar) (IdxBSRTheta idx) <br>_Convert a 2D (r,theta) bspline index into a polar bspline index._  |


























## Detailed Description


Assemble the finite element stiffness matrix with entries \(\int \alpha \nabla B_i \cdot \nabla B_j + \beta B_i B_j dx\) for \(0 < i,j < N\) where \(N\) is the number of polar splines. We only consider splines that vanish on the boundary.


For more details see the `PolarSplineFEMPoissonLikeSolver` and Emily Bourne's thesis ("Non-Uniform Numerical Schemes for the Modelling of Turbulence
in the 5D GYSELA Code". December 2022.) 

**Template parameters:**


* [**GridR**](structGridR.md) The radial grid type. 
* [**GridR**](structGridR.md) The poloidal grid type. 
* [**PolarBSplinesRTheta**](structPolarBSplinesRTheta.md) The type of the 2D polar B-splines (on the coordinate system \((r,\theta)\) including B-splines which traverse the O point). 
* `QDimRMesh` The radial quadrature grid type. 
* `QDimThetaMesh` The poloidal quadrature grid type. 
* `MetricTensorType` The type of the metric tensor associated with the mapping from the logical domain to the physical domain where the equation is defined. 
* `IdxRangeFull` The full index range of \(\phi\) including any batch dimensions. One matrix is assembled for each index in the batch dimensions. 




    
## Public Types Documentation




### typedef R 

_The radial dimension._ 
```C++
using PolarSplineFEMPoissonLikeAssembler< GridR, GridTheta, PolarBSplinesRTheta, QDimRMesh, QDimThetaMesh, MetricTensorType, IdxRangeFull >::R =  typename GridR::continuous_dimension_type;
```




<hr>



### typedef TensorBSplineStash 

_The values and derivatives of the 1D B-splines at the quadrature points._ 
```C++
using PolarSplineFEMPoissonLikeAssembler< GridR, GridTheta, PolarBSplinesRTheta, QDimRMesh, QDimThetaMesh, MetricTensorType, IdxRangeFull >::TensorBSplineStash =  detail_poisson::
            TensorBSplineStash<QDimRMesh, QDimThetaMesh, BSplinesR, BSplinesTheta>;
```




<hr>



### typedef Theta 

_The poloidal dimension._ 
```C++
using PolarSplineFEMPoissonLikeAssembler< GridR, GridTheta, PolarBSplinesRTheta, QDimRMesh, QDimThetaMesh, MetricTensorType, IdxRangeFull >::Theta =  typename GridTheta::continuous_dimension_type;
```




<hr>
## Public Functions Documentation




### function PolarSplineFEMPoissonLikeAssembler 

_Instantiate the assembler operator._ 
```C++
inline explicit PolarSplineFEMPoissonLikeAssembler::PolarSplineFEMPoissonLikeAssembler (
    Field< double, IdxRangeQuadratureRTheta > int_volume,
    MetricTensorType const & metric_tensor,
    IdxRangeFull idx_range_full
) 
```





**Parameters:**


* `int_volume` The initialised field of Jacobian values of the mapping. 
* `metric_tensor` The metric tensor associated with the mapping from the logical domain to the physical domain where the equation is defined. 
* `idx_range_full` The full index range of \(\phi\) including any batch dimensions. 




        

<hr>



### function compute\_singular\_singular\_col\_idx 

_Computes the column indices of the matrix elements corresponding to the inner products of singular basis functions and singular basis functions._ 
```C++
inline void PolarSplineFEMPoissonLikeAssembler::compute_singular_singular_col_idx (
    Kokkos::View< int *, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace > const col_idx_csr,
    Kokkos::View< int *, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace > const nnz_per_row_csr
) 
```





**Parameters:**


* `col_idx_csr` A 1D Kokkos view which stores the column indices for each non-zero component.(only for one matrix). 
* `nnz_per_row_csr` A 1D Kokkos view of length matrix\_size+1 which stores the count of the non-zeros along the lines of the matrix. 




        

<hr>



### function compute\_singular\_singular\_elements 

_Computes the matrix elements corresponding to the inner products of singular basis functions and singular basis functions._ 
```C++
template<class CoeffAlpha, class CoeffBeta>
inline void PolarSplineFEMPoissonLikeAssembler::compute_singular_singular_elements (
    CoeffAlpha const & coeff_alpha,
    CoeffBeta const & coeff_beta,
    MetricTensorType const & metric_tensor,
    Kokkos::View< double **, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace > const values_csr,
    Kokkos::View< int *, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace > const col_idx_csr,
    Kokkos::View< int *, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace > const nnz_per_row_csr
) 
```





**Parameters:**


* `coeff_alpha` A callable with signature `double operator()(IdxBatch, CoordRTheta)` returning \(\alpha\) at the given batch and coordinate. 
* `coeff_beta` A callable with signature `double operator()(IdxBatch, CoordRTheta)` returning \(\beta\) at the given batch and coordinate. 
* `metric_tensor` The metric tensor associated with the mapping from the logical domain to the physical domain where the equation is defined. 
* `values_csr` A 2D Kokkos view which stores the values of non-zero elements for the whole batch. 
* `col_idx_csr` A 1D Kokkos view which stores the column indices for each non-zero component.(only for one matrix). 
* `nnz_per_row_csr` A 1D Kokkos view of length matrix\_size+1 which stores the count of the non-zeros along the lines of the matrix. 




        

<hr>



### function compute\_singular\_tensor\_col\_idx 

_Computes the column indices of the matrix elements corresponding to the inner products of singular basis functions and tensor basis functions._ 
```C++
inline void PolarSplineFEMPoissonLikeAssembler::compute_singular_tensor_col_idx (
    Kokkos::View< int *, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace > const col_idx_csr,
    Kokkos::View< int *, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace > const nnz_per_row_csr
) 
```





**Parameters:**


* `col_idx_csr` A 1D Kokkos view which stores the column indices for each non-zero component.(only for one matrix) 
* `nnz_per_row_csr` A 1D Kokkos view of length matrix\_size+1 which stores the count of the non-zeros along the lines of the matrix. 




        

<hr>



### function compute\_singular\_tensor\_elements 

_Computes the matrix elements corresponding to the inner products of singular basis functions and tensor basis functions._ 
```C++
template<class CoeffAlpha, class CoeffBeta>
inline void PolarSplineFEMPoissonLikeAssembler::compute_singular_tensor_elements (
    CoeffAlpha const & coeff_alpha,
    CoeffBeta const & coeff_beta,
    MetricTensorType const & metric_tensor,
    Kokkos::View< double **, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace > const values_csr,
    Kokkos::View< int *, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace > const col_idx_csr,
    Kokkos::View< int *, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace > const nnz_per_row_csr
) 
```





**Parameters:**


* `coeff_alpha` A callable with signature `double operator()(IdxBatch, CoordRTheta)` returning \(\alpha\) at the given batch and coordinate. 
* `coeff_beta` A callable with signature `double operator()(IdxBatch, CoordRTheta)` returning \(\beta\) at the given batch and coordinate. 
* `metric_tensor` The metric tensor associated with the mapping from the logical domain to the physical domain where the equation is defined. 
* `values_csr` A 2D Kokkos view which stores the values of non-zero elements for the whole batch. 
* `col_idx_csr` A 1D Kokkos view which stores the column indices for each non-zero component.(only for one matrix) 
* `nnz_per_row_csr` A 1D Kokkos view of length matrix\_size+1 which stores the count of the non-zeros along the lines of the matrix. 




        

<hr>



### function compute\_tensor\_tensor\_col\_idx 

_Computes the column indices of the matrix elements corresponding to the inner products of tensor basis functions and tensor basis functions._ 
```C++
inline void PolarSplineFEMPoissonLikeAssembler::compute_tensor_tensor_col_idx (
    Kokkos::View< int *, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace > const col_idx_csr,
    Kokkos::View< int *, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace > const nnz_per_row_csr
) 
```





**Parameters:**


* `col_idx_csr` A 1D Kokkos view which stores the column indices for each non-zero component.(only for one matrix) 
* `nnz_per_row_csr` A 1D Kokkos view of length matrix\_size+1 which stores the count of the non-zeros along the lines of the matrix. 




        

<hr>



### function compute\_tensor\_tensor\_elements 

_Computes the matrix elements corresponding to the inner products of tensor basis functions and tensor basis functions._ 
```C++
template<class CoeffAlpha, class CoeffBeta>
inline void PolarSplineFEMPoissonLikeAssembler::compute_tensor_tensor_elements (
    CoeffAlpha const & coeff_alpha,
    CoeffBeta const & coeff_beta,
    MetricTensorType const & metric_tensor,
    Kokkos::View< double **, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace > const values_csr,
    Kokkos::View< int *, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace > const col_idx_csr,
    Kokkos::View< int *, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace > const nnz_per_row_csr
) 
```





**Parameters:**


* `coeff_alpha` A callable with signature `double operator()(IdxBatch, CoordRTheta)` returning \(\alpha\) at the given batch and coordinate. 
* `coeff_beta` A callable with signature `double operator()(IdxBatch, CoordRTheta)` returning \(\beta\) at the given batch and coordinate. 
* `metric_tensor` The metric tensor associated with the mapping from the logical domain to the physical domain where the equation is defined. 
* `values_csr` A 2D Kokkos view which stores the values of non-zero elements for the whole batch. 
* `col_idx_csr` A 1D Kokkos view which stores the column indices for each non-zero component.(only for one matrix) 
* `nnz_per_row_csr` A 1D Kokkos view of length matrix\_size+1 which stores the count of the non-zeros along the lines of the matrix. 




        

<hr>



### function get\_tensor\_bspline\_stash 

_Get the_ [_**TensorBSplineStash**_](classPolarSplineFEMPoissonLikeAssembler.md#typedef-tensorbsplinestash) _of the values and derivatives of the 1D B-splines at the quadrature points._
```C++
inline TensorBSplineStash PolarSplineFEMPoissonLikeAssembler::get_tensor_bspline_stash () const
```





**Returns:**

The stashed values of the 1D B-splines. 





        

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
template<class CoeffAlpha, class CoeffBeta>
inline void PolarSplineFEMPoissonLikeAssembler::operator() (
    std::unique_ptr< MatrixBatchCsr < Kokkos::DefaultExecutionSpace, MatrixBatchCsrSolver::CG > > const & gko_matrix,
    CoeffAlpha const & coeff_alpha,
    CoeffBeta const & coeff_beta
) 
```





**Parameters:**


* `gko_matrix` The pointer to the assembled matrix. 
* `coeff_alpha` A callable object with signature `double operator()(IdxBatch, CoordRTheta)` returning the value of \(\alpha\) at the given batch and coordinate. 
* `coeff_beta` A callable object with signature `double operator()(IdxBatch, CoordRTheta)` returning the value of \(\beta\) at the given batch and coordinate.



**Template parameters:**


* `CoeffAlpha` A callable type for evaluating \(\alpha\) at a batch and coordinate. 
* `CoeffBeta` A callable type for evaluating \(\beta\) at a batch and coordinate. 




        

<hr>



### function setup\_sparse\_matrix 

_Compute the sparsity pattern for the stiffness matrix._ 
```C++
inline void PolarSplineFEMPoissonLikeAssembler::setup_sparse_matrix (
    std::unique_ptr< MatrixBatchCsr < Kokkos::DefaultExecutionSpace, MatrixBatchCsrSolver::CG > > & gko_matrix,
    std::optional< int > max_iter=std::nullopt,
    std::optional< double > res_tol=std::nullopt,
    std::optional< bool > batch_solver_logger=std::nullopt,
    std::optional< int > preconditioner_max_block_size=std::nullopt,
    std::optional< GkoMatrixResidualNorm > residual_norm=std::nullopt
) 
```





**Parameters:**


* `gko_matrix` The pointer to the assembled matrix. 
* `max_iter` The maximum number of iterations possible for the batched CSR solver. 
* `res_tol` The residual tolerance for the batched CSR solver. The solver stops when \|\|r\|\|/\|\|b\|\| &lt;= res\_tol where r is the residual computed by the CG recurrence. 
* `batch_solver_logger` Indicates whether log information such as the residual and the number of iterations should be monitored. 
* `preconditioner_max_block_size` The maximum size of the Jacobi preconditioner used by the batched CSR solver. 
* `residual_norm` The norm used to evaluate the residual in the stopping criterion.



**Template parameters:**


* `Mapping` A class describing a mapping from curvilinear coordinates to Cartesian coordinates. 




        

<hr>
## Public Static Functions Documentation




### function fill\_test\_weights 

_Computes the test\_mass\_weight and test\_stiffness\_weight at every quadrature point in the support of a test function and saves them in team scratch memory._ 
```C++
template<class CoeffAlpha, class CoeffBeta>
static inline KOKKOS_FUNCTION void PolarSplineFEMPoissonLikeAssembler::fill_test_weights (
    const Kokkos::TeamPolicy<>::member_type & team,
    QuadField test_mass,
    DVectorField< IdxRangeQuadratureRTheta, VectorIndexSet< R R , Theta Theta > > test_stiffness,
    IdxBSRTheta idx_test,
    IdxBatch idx_batch,
    CoeffAlpha const & coeff_alpha,
    CoeffBeta const & coeff_beta,
    MetricTensorType const & metric_tensor,
    IdxRangeQuadratureRTheta const & full_quad_idx_range,
    DField< IdxRangeQuadratureRTheta > int_volume,
    TensorBSplineStash const & bspline_stash
) 
```





**Parameters:**


* `team` The team of threads from which this function is called. 
* `test_mass` The field where test\_mass\_weight is saved. Its index range is the range of quadrature points in the support of the test B-spline (see detail\_poisson::get\_support\_quadrature). 
* `test_stiffness` The field where test\_stiffness\_weight is saved. It is defined on the same index range as test\_mass. 
* `idx_test` The 2D index of the test B-spline. 
* `idx_batch` The index of the batch for which the element is being calculated. 
* `coeff_alpha` A callable with signature `double operator()(CoordRTheta)` returning \(\alpha\). 
* `coeff_beta` A callable with signature `double operator()(CoordRTheta)` returning \(\beta\). 
* `metric_tensor` The metric tensor associated with the mapping from the logical domain to the physical domain where the equation is defined. 
* `full_quad_idx_range` The index range of all the quadrature points. 
* `int_volume` The integral volume associated with each point used in the quadrature scheme. 
* `bspline_stash` The values and derivatives of the 1D B-splines at the quadrature points. 




        

<hr>



### function get\_matrix\_stencil\_element 

_Computes the matrix element corresponding to two tensor product splines with index idx\_test and idx\_trial._ 
```C++
static inline KOKKOS_FUNCTION double PolarSplineFEMPoissonLikeAssembler::get_matrix_stencil_element (
    const Kokkos::TeamPolicy<>::member_type & team,
    IdxBSRTheta idx_test,
    IdxBSRTheta idx_trial,
    QuadField test_mass,
    DVectorField< IdxRangeQuadratureRTheta, VectorIndexSet< R R , Theta Theta > > test_stiffness,
    IdxRangeQuadratureRTheta const & full_quad_idx_range,
    TensorBSplineStash const & bspline_stash
) 
```





**Parameters:**


* `team` The team of threads from which this function is called. 
* `idx_test` The 2D index of the tensor product B-spline in the test space. 
* `idx_trial` The 2D index of the tensor product B-spline in the trial space. 
* `test_mass` The test\_mass\_weight at the quadrature points in the support of the test B-spline (see fill\_test\_weights). 
* `test_stiffness` The test\_stiffness\_weight at the quadrature points in the support of the test B-spline (see fill\_test\_weights). 
* `full_quad_idx_range` The index range of all the quadrature points. 
* `bspline_stash` The values and derivatives of the 1D B-splines at the quadrature points. 



**Returns:**

The value of the matrix element. 





        

<hr>



### function get\_test\_weights 

_Computes the test\_stiffness\_weight and test\_mass\_weight for a given quadrature point._ 
```C++
template<class CoeffAlpha, class CoeffBeta, class IdxTestType>
static inline KOKKOS_FUNCTION std::tuple< double, DVector< R R , Theta Theta > > PolarSplineFEMPoissonLikeAssembler::get_test_weights (
    IdxTestType idx_test,
    IdxQuadratureRTheta idx_quad,
    IdxBatch idx_batch,
    CoeffAlpha const & coeff_alpha,
    CoeffBeta const & coeff_beta,
    MetricTensorType const & metric_tensor,
    DField< IdxRangeQuadratureRTheta > int_volume,
    TensorBSplineStash const & bspline_stash
) 
```



The weak integral element is \(w (\alpha \nabla B_i \cdot G^{-1} \nabla B_j + \beta B_i B_j)\) where \(w\) is the integral volume, \(B_i\) is the test function and \(B_j\) is the trial function. This function computes \(w \alpha G^{-1} \nabla B_i\) and \(w \beta B_i\) so that the weak integral element can be computed from the values and derivatives of the trial function alone.




**Parameters:**


* `idx_test` The index of the test basis spline. 
* `idx_quad` The index for the point in the quadrature scheme. 
* `idx_batch` The index of the batch for which the element is being calculated. 
* `coeff_alpha` A callable with signature `double operator()(IdxBatch, CoordRTheta)` returning \(\alpha\) at the given batch and coordinate. 
* `coeff_beta` A callable with signature `double operator()(IdxBatch, CoordRTheta)` returning \(\beta\) at the given batch and coordinate. 
* `metric_tensor` The metric tensor associated with the mapping from the logical domain to the physical domain where the equation is defined. 
* `int_volume` The integral volume associated with each point used in the quadrature scheme. 
* `bspline_stash` The values and derivatives of the 1D B-splines at the quadrature points. 



**Returns:**

A tuple containing
* test\_mass\_weight, the quantity \(w \beta B_i\) which multiplies the value of the trial function,
* test\_stiffness\_weight, the quantity \(w \alpha G^{-1} \nabla B_i\) which multiplies the gradient of the trial function. 







        

<hr>



### function stash\_basis\_eval 

_Fill the stash with the values and derivatives of the non-zero 1D B-splines at the 1D quadrature points._ 
```C++
template<class QDimMesh, class BSplines, class DerivDim>
static inline void PolarSplineFEMPoissonLikeAssembler::stash_basis_eval (
    DField< IdxRange< QDimMesh, BSplines, DerivDim > > basis_vals,
    Field< Idx< BSplines >, IdxRange< QDimMesh > > first_basis
) 
```





**Parameters:**


* `basis_vals` The values and derivatives of the non-zero 1D B-splines. The B-spline index is the offset from the first non-zero B-spline. 
* `first_basis` The first non-zero B-spline at each quadrature point.



**Template parameters:**


* `QDimMesh` The 1D quadrature grid. 
* `BSplines` The 1D B-splines. 
* `DerivDim` The derivative dimension. 




        

<hr>



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

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/pde_solvers/polar_spline_fem_poisson_like_assembler.hpp`

