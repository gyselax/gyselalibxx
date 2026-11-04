

# Class TensorProductBSplinePoissonAssembler

**template &lt;typename Grid1, typename Grid2, typename BSplines1, typename BSplines2, typename QDim1Mesh, typename QDim2Mesh, concepts::MetricTensor MetricTensorType, typename IdxRangeFull&gt;**



[**ClassList**](annotated.md) **>** [**TensorProductBSplinePoissonAssembler**](classTensorProductBSplinePoissonAssembler.md)



_An operator to assemble a Poisson-like stiffness matrix using tensor product B-splines._ [More...](#detailed-description)

* `#include <tensor_product_bspline_poisson_assembler.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef typename Grid1::continuous\_dimension\_type | [**Dim1**](#typedef-dim1)  <br>_The first dimension._  |
| typedef typename Grid2::continuous\_dimension\_type | [**Dim2**](#typedef-dim2)  <br>_The second dimension._  |
| typedef IdxRange&lt; BSplines1, BSplines2 &gt; | [**IdxRangeBS2D**](#typedef-idxrangebs2d)  <br>_The index range of the tensor product B-splines._  |
| typedef ddc::remove\_dims\_of\_t&lt; IdxRangeFull, Grid1, Grid2 &gt; | [**IdxRangeBatch**](#typedef-idxrangebatch)  <br>_The index range of the batch dimensions._  |
| typedef IdxRange&lt; QDim1Mesh, QDim2Mesh &gt; | [**IdxRangeQuadrature2D**](#typedef-idxrangequadrature2d)  <br>_The index range of the quadrature points._  |
| typedef [**detail\_poisson::TensorBSplineStash**](structdetail__poisson_1_1TensorBSplineStash.md)&lt; QDim1Mesh, QDim2Mesh, BSplines1, BSplines2 &gt; | [**TensorBSplineStash**](#typedef-tensorbsplinestash)  <br>_The values and derivatives of the 1D B-splines at the quadrature points._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**TensorProductBSplinePoissonAssembler**](#function-tensorproductbsplinepoissonassembler) (QuadField int\_volume, MetricTensorType const & metric\_tensor, IdxRangeFull idx\_range\_full, [**IdxRangeBS2D**](classTensorProductBSplinePoissonAssembler.md#typedef-idxrangebs2d) idx\_range\_tensor\_basis) <br>_Instantiate the assembler operator._  |
|  void | [**compute\_col\_idx**](#function-compute_col_idx) (Kokkos::View&lt; int \*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace &gt; const col\_idx\_csr, Kokkos::View&lt; int \*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace &gt; const nnz\_per\_row\_csr, int first\_col\_offset) const<br>_Computes the column indices of the matrix elements corresponding to the inner products of tensor basis functions and tensor basis functions._  |
|  void | [**compute\_elements**](#function-compute_elements) (CoeffAlpha const & coeff\_alpha, CoeffBeta const & coeff\_beta, Kokkos::View&lt; double \*\*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace &gt; const values\_csr, Kokkos::View&lt; int \*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace &gt; const col\_idx\_csr, Kokkos::View&lt; int \*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace &gt; const nnz\_per\_row\_csr, int first\_col\_offset) const<br>_Computes the matrix elements corresponding to the inner products of tensor basis functions and tensor basis functions._  |
|  QuadField | [**get\_int\_volume**](#function-get_int_volume) () const<br>_Get the integral volume associated with each quadrature point._  |
|  [**IdxRangeBS2D**](classTensorProductBSplinePoissonAssembler.md#typedef-idxrangebs2d) | [**get\_tensor\_basis\_idx\_range**](#function-get_tensor_basis_idx_range) () const<br>_Get the index range of the tensor product B-splines in the FEM basis._  |
|  [**TensorBSplineStash**](classTensorProductBSplinePoissonAssembler.md#typedef-tensorbsplinestash) | [**get\_tensor\_bspline\_stash**](#function-get_tensor_bspline_stash) () const<br>_Get the_ [_**TensorBSplineStash**_](classTensorProductBSplinePoissonAssembler.md#typedef-tensorbsplinestash) _of the values and derivatives of the 1D B-splines at the quadrature points._ |
|  void | [**init\_nnz\_per\_line**](#function-init_nnz_per_line) (Kokkos::View&lt; int \*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace &gt; nnz\_per\_row) const<br>_Add the number of non-zero elements per line of the tensor-tensor matrix to a nnz data structure._  |
|  void | [**init\_tensor\_tensor\_nnz\_per\_line**](#function-init_tensor_tensor_nnz_per_line) (Kokkos::View&lt; int \*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace &gt; nnz\_per\_row) const<br>_Fills the nnz data structure of the matrix corresponding to the inner products of tensor basis functions and tensor basis functions._  |
|  int | [**n\_tensor\_tensor\_elements**](#function-n_tensor_tensor_elements) () const<br>_Get the number of non-zero elements in the matrix corresponding to the inner products of tensor basis functions and tensor basis functions._  |
|  void | [**operator()**](#function-operator) (std::unique\_ptr&lt; [**MatrixBatchCsr**](classMatrixBatchCsr.md)&lt; Kokkos::DefaultExecutionSpace, MatrixBatchCsrSolver::CG &gt; &gt; const & gko\_matrix, CoeffAlpha const & coeff\_alpha, CoeffBeta const & coeff\_beta) const<br>_Assemble the stiffness matrix._  |
|  void | [**setup\_sparse\_matrix**](#function-setup_sparse_matrix) (std::unique\_ptr&lt; [**MatrixBatchCsr**](classMatrixBatchCsr.md)&lt; Kokkos::DefaultExecutionSpace, MatrixBatchCsrSolver::CG &gt; &gt; & gko\_matrix, std::optional&lt; int &gt; max\_iter=std::nullopt, std::optional&lt; double &gt; res\_tol=std::nullopt, std::optional&lt; bool &gt; batch\_solver\_logger=std::nullopt, std::optional&lt; int &gt; preconditioner\_max\_block\_size=std::nullopt, std::optional&lt; GkoMatrixResidualNorm &gt; residual\_norm=std::nullopt) const<br>_Compute the sparsity pattern for the stiffness matrix._  |


## Public Static Functions

| Type | Name |
| ---: | :--- |
|  KOKKOS\_FUNCTION void | [**fill\_test\_weights**](#function-fill_test_weights) (const Kokkos::TeamPolicy&lt;&gt;::member\_type & team, QuadField test\_mass, DVectorField&lt; [**IdxRangeQuadrature2D**](classTensorProductBSplinePoissonAssembler.md#typedef-idxrangequadrature2d), VectorIndexSet&lt; [**Dim1**](classTensorProductBSplinePoissonAssembler.md#typedef-dim1), [**Dim2**](classTensorProductBSplinePoissonAssembler.md#typedef-dim2) &gt; &gt; test\_stiffness, IdxBS2D idx\_test, IdxBatch idx\_batch, CoeffAlpha const & coeff\_alpha, CoeffBeta const & coeff\_beta, MetricTensorType const & metric\_tensor, [**IdxRangeQuadrature2D**](classTensorProductBSplinePoissonAssembler.md#typedef-idxrangequadrature2d) const & full\_quad\_idx\_range, DField&lt; [**IdxRangeQuadrature2D**](classTensorProductBSplinePoissonAssembler.md#typedef-idxrangequadrature2d) &gt; int\_volume, [**TensorBSplineStash**](classTensorProductBSplinePoissonAssembler.md#typedef-tensorbsplinestash) const & bspline\_stash) <br>_Computes the test\_mass\_weight and test\_stiffness\_weight at every quadrature point in the support of a test function and saves them in team scratch memory._  |
|  KOKKOS\_FUNCTION double | [**get\_matrix\_stencil\_element**](#function-get_matrix_stencil_element) (const Kokkos::TeamPolicy&lt;&gt;::member\_type & team, IdxBS2D idx\_test, IdxBS2D idx\_trial, QuadField test\_mass, DVectorField&lt; [**IdxRangeQuadrature2D**](classTensorProductBSplinePoissonAssembler.md#typedef-idxrangequadrature2d), VectorIndexSet&lt; [**Dim1**](classTensorProductBSplinePoissonAssembler.md#typedef-dim1), [**Dim2**](classTensorProductBSplinePoissonAssembler.md#typedef-dim2) &gt; &gt; test\_stiffness, [**IdxRangeQuadrature2D**](classTensorProductBSplinePoissonAssembler.md#typedef-idxrangequadrature2d) const & full\_quad\_idx\_range, [**TensorBSplineStash**](classTensorProductBSplinePoissonAssembler.md#typedef-tensorbsplinestash) const & bspline\_stash) <br>_Computes the matrix element corresponding to two tensor product splines with index idx\_test and idx\_trial._  |
|  KOKKOS\_FUNCTION std::tuple&lt; double, DVector&lt; [**Dim1**](classTensorProductBSplinePoissonAssembler.md#typedef-dim1), [**Dim2**](classTensorProductBSplinePoissonAssembler.md#typedef-dim2) &gt; &gt; | [**get\_test\_weights**](#function-get_test_weights) (IdxBS2D idx\_test, IdxQuadrature2D idx\_quad, IdxBatch idx\_batch, CoeffAlpha const & coeff\_alpha, CoeffBeta const & coeff\_beta, MetricTensorType const & metric\_tensor, DField&lt; [**IdxRangeQuadrature2D**](classTensorProductBSplinePoissonAssembler.md#typedef-idxrangequadrature2d) &gt; int\_volume, [**TensorBSplineStash**](classTensorProductBSplinePoissonAssembler.md#typedef-tensorbsplinestash) const & bspline\_stash) <br>_Computes the test\_stiffness\_weight and test\_mass\_weight for a given quadrature point for a tensor product test basis spline._  |
|  void | [**stash\_basis\_eval**](#function-stash_basis_eval) (DField&lt; IdxRange&lt; QDimMesh, BSplines, DerivDim &gt; &gt; basis\_vals, Field&lt; Idx&lt; BSplines &gt;, IdxRange&lt; QDimMesh &gt; &gt; first\_basis) <br>_Fill the stash with the values and derivatives of the non-zero 1D B-splines at the 1D quadrature points._  |


























## Detailed Description


The stiffness matrix has entries \(\int \alpha \nabla B_i \cdot \nabla B_j + \beta B_i B_j dx\) where \(B_i\) and \(B_j\) are tensor product B-splines.


The B-splines and the quadrature points both have a tensor product structure. The values and derivatives of the 1D B-splines at the 1D quadrature points are therefore stored (see [**detail\_poisson::TensorBSplineStash**](structdetail__poisson_1_1TensorBSplineStash.md)) so that the tensor product B-splines can be evaluated at the quadrature points without being recalculated.


Each dimension may be periodic or non-periodic (as determined by the B-splines). In a periodic dimension the basis must contain at least 2\*degree+1 B-splines so that a B-spline does not overlap with itself across the boundary.




**Template parameters:**


* `Grid1` The grid type in the first dimension. 
* `Grid2` The grid type in the second dimension. 
* `BSplines1` The 1D B-splines in the first dimension. 
* `BSplines2` The 1D B-splines in the second dimension. 
* `QDim1Mesh` The quadrature grid type in the first dimension. 
* `QDim2Mesh` The quadrature grid type in the second dimension. 
* `MetricTensorType` The type of the metric tensor associated with the mapping from the logical domain to the physical domain where the equation is defined. 
* `IdxRangeFull` The full index range of \(\phi\) including any batch dimensions. One matrix is assembled for each index in the batch dimensions. 




    
## Public Types Documentation




### typedef Dim1 

_The first dimension._ 
```C++
using TensorProductBSplinePoissonAssembler< Grid1, Grid2, BSplines1, BSplines2, QDim1Mesh, QDim2Mesh, MetricTensorType, IdxRangeFull >::Dim1 =  typename Grid1::continuous_dimension_type;
```




<hr>



### typedef Dim2 

_The second dimension._ 
```C++
using TensorProductBSplinePoissonAssembler< Grid1, Grid2, BSplines1, BSplines2, QDim1Mesh, QDim2Mesh, MetricTensorType, IdxRangeFull >::Dim2 =  typename Grid2::continuous_dimension_type;
```




<hr>



### typedef IdxRangeBS2D 

_The index range of the tensor product B-splines._ 
```C++
using TensorProductBSplinePoissonAssembler< Grid1, Grid2, BSplines1, BSplines2, QDim1Mesh, QDim2Mesh, MetricTensorType, IdxRangeFull >::IdxRangeBS2D =  IdxRange<BSplines1, BSplines2>;
```




<hr>



### typedef IdxRangeBatch 

_The index range of the batch dimensions._ 
```C++
using TensorProductBSplinePoissonAssembler< Grid1, Grid2, BSplines1, BSplines2, QDim1Mesh, QDim2Mesh, MetricTensorType, IdxRangeFull >::IdxRangeBatch =  ddc::remove_dims_of_t<IdxRangeFull, Grid1, Grid2>;
```




<hr>



### typedef IdxRangeQuadrature2D 

_The index range of the quadrature points._ 
```C++
using TensorProductBSplinePoissonAssembler< Grid1, Grid2, BSplines1, BSplines2, QDim1Mesh, QDim2Mesh, MetricTensorType, IdxRangeFull >::IdxRangeQuadrature2D =  IdxRange<QDim1Mesh, QDim2Mesh>;
```




<hr>



### typedef TensorBSplineStash 

_The values and derivatives of the 1D B-splines at the quadrature points._ 
```C++
using TensorProductBSplinePoissonAssembler< Grid1, Grid2, BSplines1, BSplines2, QDim1Mesh, QDim2Mesh, MetricTensorType, IdxRangeFull >::TensorBSplineStash =  detail_poisson::TensorBSplineStash<QDim1Mesh, QDim2Mesh, BSplines1, BSplines2>;
```




<hr>
## Public Functions Documentation




### function TensorProductBSplinePoissonAssembler 

_Instantiate the assembler operator._ 
```C++
inline explicit TensorProductBSplinePoissonAssembler::TensorProductBSplinePoissonAssembler (
    QuadField int_volume,
    MetricTensorType const & metric_tensor,
    IdxRangeFull idx_range_full,
    IdxRangeBS2D idx_range_tensor_basis
) 
```





**Parameters:**


* `int_volume` The initialised field of Jacobian values of the mapping. 
* `metric_tensor` The metric tensor associated with the mapping from the logical domain to the physical domain where the equation is defined. 
* `idx_range_full` The full index range of \(\phi\) including any batch dimensions. 
* `idx_range_tensor_basis` The index range of the tensor product B-splines in the FEM basis. B-splines which are excluded from this index range (e.g. to impose boundary conditions) are not included in the matrix. In a periodic dimension all B-splines must be included. 




        

<hr>



### function compute\_col\_idx 

_Computes the column indices of the matrix elements corresponding to the inner products of tensor basis functions and tensor basis functions._ 
```C++
inline void TensorProductBSplinePoissonAssembler::compute_col_idx (
    Kokkos::View< int *, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace > const col_idx_csr,
    Kokkos::View< int *, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace > const nnz_per_row_csr,
    int first_col_offset
) const
```



The tensor-tensor elements are found at the end of each row of the tensor-tensor matrix.




**Parameters:**


* `col_idx_csr` A 1D Kokkos view which stores the column indices for each non-zero component of the full matrix (only for one matrix). 
* `nnz_per_row_csr` A 1D Kokkos view of length n+1, where n is the number of tensor basis functions, containing the elements of the nnz data structure of the full matrix which correspond to the rows of the tensor-tensor matrix. 
* `first_col_offset` The index of the first column of the tensor-tensor matrix in the full matrix. 




        

<hr>



### function compute\_elements 

_Computes the matrix elements corresponding to the inner products of tensor basis functions and tensor basis functions._ 
```C++
template<class CoeffAlpha, class CoeffBeta>
inline void TensorProductBSplinePoissonAssembler::compute_elements (
    CoeffAlpha const & coeff_alpha,
    CoeffBeta const & coeff_beta,
    Kokkos::View< double **, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace > const values_csr,
    Kokkos::View< int *, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace > const col_idx_csr,
    Kokkos::View< int *, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace > const nnz_per_row_csr,
    int first_col_offset
) const
```





**Parameters:**


* `coeff_alpha` A callable with signature `double operator()(IdxBatch, Coord< Dim1 , Dim2 >)` returning \(\alpha\) at the given batch and coordinate. 
* `coeff_beta` A callable with signature `double operator()(IdxBatch, Coord< Dim1 , Dim2 >)` returning \(\beta\) at the given batch and coordinate. 
* `values_csr` A 2D Kokkos view which stores the values of non-zero elements of the full matrix for the whole batch. 
* `col_idx_csr` A 1D Kokkos view which stores the column indices for each non-zero component of the full matrix (only for one matrix). 
* `nnz_per_row_csr` A 1D Kokkos view of length n+1, where n is the number of tensor basis functions, containing the elements of the nnz data structure of the full matrix which correspond to the rows of the tensor-tensor matrix. 
* `first_col_offset` The index of the first column of the tensor-tensor matrix in the full matrix. 




        

<hr>



### function get\_int\_volume 

_Get the integral volume associated with each quadrature point._ 
```C++
inline QuadField TensorProductBSplinePoissonAssembler::get_int_volume () const
```





**Returns:**

The integral volume. 





        

<hr>



### function get\_tensor\_basis\_idx\_range 

_Get the index range of the tensor product B-splines in the FEM basis._ 
```C++
inline IdxRangeBS2D TensorProductBSplinePoissonAssembler::get_tensor_basis_idx_range () const
```





**Returns:**

The index range of the tensor product B-splines. 





        

<hr>



### function get\_tensor\_bspline\_stash 

_Get the_ [_**TensorBSplineStash**_](classTensorProductBSplinePoissonAssembler.md#typedef-tensorbsplinestash) _of the values and derivatives of the 1D B-splines at the quadrature points._
```C++
inline TensorBSplineStash TensorProductBSplinePoissonAssembler::get_tensor_bspline_stash () const
```





**Returns:**

The stashed values of the 1D B-splines. 





        

<hr>



### function init\_nnz\_per\_line 

_Add the number of non-zero elements per line of the tensor-tensor matrix to a nnz data structure._ 
```C++
inline void TensorProductBSplinePoissonAssembler::init_nnz_per_line (
    Kokkos::View< int *, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace > nnz_per_row
) const
```



After this function the array will contain: nnz\_per\_row[0] = unchanged (0 for a standalone matrix). nnz\_per\_row[1] += number of non-zero elements in line 0 of the tensor-tensor matrix. ... nnz\_per\_row[n] += number of non-zero elements in lines 0-(n-1) of the tensor-tensor matrix.




**Parameters:**


* `nnz_per_row` A 1D Kokkos view of length n+1, where n is the number of tensor basis functions, containing the elements of the nnz data structure which correspond to the rows of the tensor-tensor matrix. 




        

<hr>



### function init\_tensor\_tensor\_nnz\_per\_line 

_Fills the nnz data structure of the matrix corresponding to the inner products of tensor basis functions and tensor basis functions._ 
```C++
inline void TensorProductBSplinePoissonAssembler::init_tensor_tensor_nnz_per_line (
    Kokkos::View< int *, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace > nnz_per_row
) const
```



The rows and columns of the matrix are the tensor basis functions, ordered as in the 2D index range of the tensor basis functions. After this function the array will contain: nnz\_per\_row[0] = unchanged (0 for a standalone matrix). nnz\_per\_row[1] = number of non-zero elements in line 0 of the matrix. ... nnz\_per\_row[n] = number of non-zero elements in lines 0-(n-1) of the matrix. The values are added to the values already saved in nnz\_per\_row.




**Parameters:**


* `nnz_per_row` A 1D Kokkos view of length n+1, where n is the number of tensor basis functions, which stores the sum of the non-zeros in the matrix on all lines up to the one in. 




        

<hr>



### function n\_tensor\_tensor\_elements 

_Get the number of non-zero elements in the matrix corresponding to the inner products of tensor basis functions and tensor basis functions._ 
```C++
inline int TensorProductBSplinePoissonAssembler::n_tensor_tensor_elements () const
```





**Returns:**

The number of non-zero elements in the tensor-tensor matrix. 





        

<hr>



### function operator() 

_Assemble the stiffness matrix._ 
```C++
template<class CoeffAlpha, class CoeffBeta>
inline void TensorProductBSplinePoissonAssembler::operator() (
    std::unique_ptr< MatrixBatchCsr < Kokkos::DefaultExecutionSpace, MatrixBatchCsrSolver::CG > > const & gko_matrix,
    CoeffAlpha const & coeff_alpha,
    CoeffBeta const & coeff_beta
) const
```





**Parameters:**


* `gko_matrix` The pointer to the matrix initialised by setup\_sparse\_matrix. 
* `coeff_alpha` A callable object with signature `double operator()(IdxBatch, Coord< Dim1 , Dim2 >)` returning the value of \(\alpha\) at the given batch and coordinate. 
* `coeff_beta` A callable object with signature `double operator()(IdxBatch, Coord< Dim1 , Dim2 >)` returning the value of \(\beta\) at the given batch and coordinate.



**Template parameters:**


* `CoeffAlpha` A callable type for evaluating \(\alpha\) at a batch and coordinate. 
* `CoeffBeta` A callable type for evaluating \(\beta\) at a batch and coordinate. 




        

<hr>



### function setup\_sparse\_matrix 

_Compute the sparsity pattern for the stiffness matrix._ 
```C++
inline void TensorProductBSplinePoissonAssembler::setup_sparse_matrix (
    std::unique_ptr< MatrixBatchCsr < Kokkos::DefaultExecutionSpace, MatrixBatchCsrSolver::CG > > & gko_matrix,
    std::optional< int > max_iter=std::nullopt,
    std::optional< double > res_tol=std::nullopt,
    std::optional< bool > batch_solver_logger=std::nullopt,
    std::optional< int > preconditioner_max_block_size=std::nullopt,
    std::optional< GkoMatrixResidualNorm > residual_norm=std::nullopt
) const
```



The rows and columns of the matrix are the tensor product B-splines in the FEM basis ordered following the layout\_right order of their index range.




**Parameters:**


* `gko_matrix` The pointer to the assembled matrix. 
* `max_iter` The maximum number of iterations possible for the batched CSR solver. 
* `res_tol` The residual tolerance for the batched CSR solver. The solver stops when \|\|r\|\|/\|\|b\|\| &lt;= res\_tol where r is the residual computed by the CG recurrence. 
* `batch_solver_logger` Indicates whether log information such as the residual and the number of iterations should be monitored. 
* `preconditioner_max_block_size` The maximum size of the Jacobi preconditioner used by the batched CSR solver. 
* `residual_norm` The norm used to evaluate the residual in the stopping criterion. 




        

<hr>
## Public Static Functions Documentation




### function fill\_test\_weights 

_Computes the test\_mass\_weight and test\_stiffness\_weight at every quadrature point in the support of a test function and saves them in team scratch memory._ 
```C++
template<class CoeffAlpha, class CoeffBeta>
static inline KOKKOS_FUNCTION void TensorProductBSplinePoissonAssembler::fill_test_weights (
    const Kokkos::TeamPolicy<>::member_type & team,
    QuadField test_mass,
    DVectorField< IdxRangeQuadrature2D , VectorIndexSet< Dim1 , Dim2 > > test_stiffness,
    IdxBS2D idx_test,
    IdxBatch idx_batch,
    CoeffAlpha const & coeff_alpha,
    CoeffBeta const & coeff_beta,
    MetricTensorType const & metric_tensor,
    IdxRangeQuadrature2D const & full_quad_idx_range,
    DField< IdxRangeQuadrature2D > int_volume,
    TensorBSplineStash const & bspline_stash
) 
```





**Parameters:**


* `team` The team of threads from which this function is called. 
* `test_mass` The field where test\_mass\_weight is saved. Its index range is the range of quadrature points in the support of the test B-spline (see detail\_poisson::get\_support\_quadrature). 
* `test_stiffness` The field where test\_stiffness\_weight is saved. It is defined on the same index range as test\_mass. 
* `idx_test` The 2D index of the test B-spline. 
* `idx_batch` The index of the batch for which the element is being calculated. 
* `coeff_alpha` A callable with signature `double operator()(IdxBatch, Coord< Dim1 , Dim2 >)` returning \(\alpha\) at the given batch and coordinate. 
* `coeff_beta` A callable with signature `double operator()(IdxBatch, Coord< Dim1 , Dim2 >)` returning \(\beta\) at the given batch and coordinate. 
* `metric_tensor` The metric tensor associated with the mapping from the logical domain to the physical domain where the equation is defined. 
* `full_quad_idx_range` The index range of all the quadrature points. 
* `int_volume` The integral volume associated with each point used in the quadrature scheme. 
* `bspline_stash` The values and derivatives of the 1D B-splines at the quadrature points. 




        

<hr>



### function get\_matrix\_stencil\_element 

_Computes the matrix element corresponding to two tensor product splines with index idx\_test and idx\_trial._ 
```C++
static inline KOKKOS_FUNCTION double TensorProductBSplinePoissonAssembler::get_matrix_stencil_element (
    const Kokkos::TeamPolicy<>::member_type & team,
    IdxBS2D idx_test,
    IdxBS2D idx_trial,
    QuadField test_mass,
    DVectorField< IdxRangeQuadrature2D , VectorIndexSet< Dim1 , Dim2 > > test_stiffness,
    IdxRangeQuadrature2D const & full_quad_idx_range,
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

_Computes the test\_stiffness\_weight and test\_mass\_weight for a given quadrature point for a tensor product test basis spline._ 
```C++
template<class CoeffAlpha, class CoeffBeta>
static inline KOKKOS_FUNCTION std::tuple< double, DVector< Dim1 , Dim2 > > TensorProductBSplinePoissonAssembler::get_test_weights (
    IdxBS2D idx_test,
    IdxQuadrature2D idx_quad,
    IdxBatch idx_batch,
    CoeffAlpha const & coeff_alpha,
    CoeffBeta const & coeff_beta,
    MetricTensorType const & metric_tensor,
    DField< IdxRangeQuadrature2D > int_volume,
    TensorBSplineStash const & bspline_stash
) 
```



See detail\_poisson::get\_test\_weights for more details.




**Parameters:**


* `idx_test` The 2D index of the test tensor product B-spline. 
* `idx_quad` The index for the point in the quadrature scheme. 
* `idx_batch` The index of the batch for which the element is being calculated. 
* `coeff_alpha` A callable with signature `double operator()(IdxBatch, Coord< Dim1 , Dim2 >)` returning \(\alpha\) at the given batch and coordinate. 
* `coeff_beta` A callable with signature `double operator()(IdxBatch, Coord< Dim1 , Dim2 >)` returning \(\beta\) at the given batch and coordinate. 
* `metric_tensor` The metric tensor associated with the mapping from the logical domain to the physical domain where the equation is defined. 
* `int_volume` The integral volume associated with each point used in the quadrature scheme. 
* `bspline_stash` The values and derivatives of the 1D B-splines at the quadrature points. 



**Returns:**

A tuple containing test\_mass\_weight and test\_stiffness\_weight. 





        

<hr>



### function stash\_basis\_eval 

_Fill the stash with the values and derivatives of the non-zero 1D B-splines at the 1D quadrature points._ 
```C++
template<class QDimMesh, class BSplines, class DerivDim>
static inline void TensorProductBSplinePoissonAssembler::stash_basis_eval (
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

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/pde_solvers/tensor_product_bspline_poisson_assembler.hpp`

