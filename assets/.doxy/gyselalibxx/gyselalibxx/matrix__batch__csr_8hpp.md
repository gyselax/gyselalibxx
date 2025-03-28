

# File matrix\_batch\_csr.hpp



[**FileList**](files.md) **>** [**matrix\_tools**](dir_8cedd1260cc2f2819c8df2fc66ad98b5.md) **>** [**matrix\_batch\_csr.hpp**](matrix__batch__csr_8hpp.md)

[Go to the source code of this file](matrix__batch__csr_8hpp_source.md)



* `#include <ginkgo/extensions/kokkos.hpp>`
* `#include <ginkgo/ginkgo.hpp>`
* `#include <Kokkos_Core.hpp>`
* `#include "matrix_batch.hpp"`
* `#include "matrix_utils.hpp"`















## Classes

| Type | Name |
| ---: | :--- |
| class | [**MatrixBatchCsr**](classMatrixBatchCsr.md) &lt;class ExecSpace, Solver&gt;<br>[_**Matrix**_](classMatrix.md) _class which is able to manage and solve a batch of sparse linear systems. Executes on either CPU or GPU. It takes advantage of the sparse structure, and the only batched solver available in Ginkgo : Stabilised Bicg. This class uses the CSR storage format which needs three arrays, one stores values, the other column indices. The third array contains the count of non-zero inside the matrix lines.(eg:for a given line index i nn\_per\_row[i]= sum of non-zeros until line i) The class returns these arrays (as Kokkos views) with the get\_csr\_views function, it is then possible to fill them outside the class. The sparsity pattern is the same for all matrices, hence column indices are stored only for one system. Tolerance and maximal number of iterations, which are parameters for the iterative solver, are set in the constructor. It is possible to get convergence information by activating the logger at constructor call._ |


## Public Types

| Type | Name |
| ---: | :--- |
| enum  | [**MatrixBatchCsrSolver**](#enum-matrixbatchcsrsolver)  <br>_A tag to choose between the batched iterative solvers provided by Ginkgo._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
|  void | [**convert\_coo\_to\_csr**](#function-convert_coo_to_csr) (std::unique\_ptr&lt; [**MatrixBatchCsr**](classMatrixBatchCsr.md)&lt; Kokkos::DefaultExecutionSpace, Solver &gt; &gt; & matrix, Kokkos::View&lt; double \*, Kokkos::LayoutRight, Kokkos::DefaultHostExecutionSpace &gt; vals\_coo\_host, Kokkos::View&lt; int \*, Kokkos::LayoutRight, Kokkos::DefaultHostExecutionSpace &gt; row\_coo\_host, Kokkos::View&lt; int \*, Kokkos::LayoutRight, Kokkos::DefaultHostExecutionSpace &gt; col\_coo\_host) <br>_A function which converts COO data storage into CSR. The results of the conversion are directly stored inside_ [_**MatrixBatchCsr**_](classMatrixBatchCsr.md) _instance and used for computations._ |




























## Public Types Documentation




### enum MatrixBatchCsrSolver 

_A tag to choose between the batched iterative solvers provided by Ginkgo._ 
```C++
enum MatrixBatchCsrSolver {
    CG,
    BICGSTAB,
    BATCH_CG,
    BATCH_BICGSTAB
};
```



BICGSTAB BiConjugate [**Gradient**](classGradient.md) Stabilise method. BICGSTAB is able to solve general sparse matrices problems. CG Conjugate [**Gradient**](classGradient.md) method (For symmetric positive definite matrices). If the matrices structures verify CG requirements, we encourage the use of this method for its lower computation cost. 


        

<hr>
## Public Functions Documentation




### function convert\_coo\_to\_csr 

_A function which converts COO data storage into CSR. The results of the conversion are directly stored inside_ [_**MatrixBatchCsr**_](classMatrixBatchCsr.md) _instance and used for computations._
```C++
template<MatrixBatchCsrSolver Solver>
void convert_coo_to_csr (
    std::unique_ptr< MatrixBatchCsr < Kokkos::DefaultExecutionSpace, Solver > > & matrix,
    Kokkos::View< double *, Kokkos::LayoutRight, Kokkos::DefaultHostExecutionSpace > vals_coo_host,
    Kokkos::View< int *, Kokkos::LayoutRight, Kokkos::DefaultHostExecutionSpace > row_coo_host,
    Kokkos::View< int *, Kokkos::LayoutRight, Kokkos::DefaultHostExecutionSpace > col_coo_host
) 
```





**Parameters:**


* `matrix` A [**MatrixBatchCsr**](classMatrixBatchCsr.md) instance, provides Kokkos views to fill. All input data is stored in 1D host allocated Kokkos View. 
* `vals_coo_host` Values. 
* `row_coo_host` Rows indices. 
* `col_coo_host` Columns indices 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/matrix_tools/matrix_batch_csr.hpp`

