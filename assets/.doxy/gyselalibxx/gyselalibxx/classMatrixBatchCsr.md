

# Class MatrixBatchCsr

**template &lt;class ExecSpace, MatrixBatchCsrSolver Solver&gt;**



[**ClassList**](annotated.md) **>** [**MatrixBatchCsr**](classMatrixBatchCsr.md)



[_**Matrix**_](classMatrix.md) _class which is able to manage and solve a batch of sparse linear systems. Executes on either CPU or GPU. It takes advantage of the sparse structure, and the only batched solver available in Ginkgo : Stabilised Bicg. This class uses the CSR storage format which needs three arrays, one stores values, the other column indices. The third array contains the count of non-zero inside the matrix lines.(eg:for a given line index i nn\_per\_row[i]= sum of non-zeros until line i) The class returns these arrays (as Kokkos views) with the get\_csr\_views function, it is then possible to fill them outside the class. The sparsity pattern is the same for all matrices, hence column indices are stored only for one system. Tolerance and maximal number of iterations, which are parameters for the iterative solver, are set in the constructor. It is possible to get convergence information by activating the logger at constructor call._[More...](#detailed-description)

* `#include <matrix_batch_csr.hpp>`



Inherits the following classes: [MatrixBatch](classMatrixBatch.md)
















## Public Types inherited from MatrixBatch

See [MatrixBatch](classMatrixBatch.md)

| Type | Name |
| ---: | :--- |
| typedef Kokkos::View&lt; double \*\*, Kokkos::LayoutRight, ExecSpace &gt; | [**BatchedRHS**](classMatrixBatch.md#typedef-batchedrhs)  <br>_The type of a Kokkos::View storing batched right-hand sides. Second dimension is batch dimension._  |






































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**MatrixBatchCsr**](#function-matrixbatchcsr-12) (const int batch\_size, const int mat\_size, const int nnz\_per\_system, std::optional&lt; int &gt; max\_iter=std::nullopt, std::optional&lt; double &gt; res\_tol=std::nullopt, std::optional&lt; bool &gt; logger=std::nullopt, std::optional&lt; int &gt; preconditionner\_max\_block\_size=1u) <br>_The constructor for_ [_**MatrixBatchCsr**_](classMatrixBatchCsr.md) _class._ |
|   | [**MatrixBatchCsr**](#function-matrixbatchcsr-22) (Kokkos::View&lt; double \*\*, Kokkos::LayoutRight, ExecSpace &gt; batch\_values, Kokkos::View&lt; int \*, Kokkos::LayoutRight, ExecSpace &gt; cols\_idx, Kokkos::View&lt; int \*, Kokkos::LayoutRight, ExecSpace &gt; nnz\_per\_row, std::optional&lt; int &gt; max\_iter=std::nullopt, std::optional&lt; double &gt; res\_tol=std::nullopt, std::optional&lt; bool &gt; logger=std::nullopt, std::optional&lt; int &gt; preconditionner\_max\_block\_size=1u) <br>_Constructor for_ [_**MatrixBatchCsr**_](classMatrixBatchCsr.md) _class._ |
|  std::tuple&lt; Kokkos::View&lt; double \*\*, Kokkos::LayoutRight, ExecSpace &gt;, Kokkos::View&lt; int \*, Kokkos::LayoutRight, ExecSpace &gt;, Kokkos::View&lt; int \*, Kokkos::LayoutRight, ExecSpace &gt; &gt; | [**get\_batch\_csr**](#function-get_batch_csr) () <br>_A function to update information about values,indices and the number of non-zero per row for the whole batch. Data is managed by Kokkos Views stored on the host._  |
|  double | [**norm**](#function-norm) (int batch\_idx) const<br>_A function returning the norm of a matrix located at batch\_idx._  |
| virtual void | [**setup\_solver**](#function-setup_solver) () <br>_Perform a pre-process operation on the solver. Must be called after filling the matrix._  |
| virtual void | [**solve**](#function-solve-12) ([**BatchedRHS**](classMatrixBatch.md#typedef-batchedrhs) const b) const<br>_Solve the batched linear problem Ax=b._  |
|  void | [**solve**](#function-solve-22) ([**BatchedRHS**](classMatrixBatch.md#typedef-batchedrhs) const x, [**BatchedRHS**](classMatrixBatch.md#typedef-batchedrhs) const b) const<br>_Solve the batched linear problem Ax=b._  |


## Public Functions inherited from MatrixBatch

See [MatrixBatch](classMatrixBatch.md)

| Type | Name |
| ---: | :--- |
|  std::size\_t | [**batch\_size**](classMatrixBatch.md#function-batch_size) () const<br>_Get the batch size of the linear problem._  |
| virtual void | [**setup\_solver**](classMatrixBatch.md#function-setup_solver) () = 0<br>_Perform a pre-process operation on the solver. Must be called after filling the matrix._  |
|  std::size\_t | [**size**](classMatrixBatch.md#function-size) () const<br>_Get the size of the square matrix corresponding to a single batch in one of its dimensions._  |
| virtual void | [**solve**](classMatrixBatch.md#function-solve) ([**BatchedRHS**](classMatrixBatch.md#typedef-batchedrhs) b) const = 0<br>_Solve the multiple right-hand sides linear problem Ax=b inplace._  |
| virtual  | [**~MatrixBatch**](classMatrixBatch.md#function-matrixbatch) () = default<br>_Destruct._  |
















































## Protected Functions inherited from MatrixBatch

See [MatrixBatch](classMatrixBatch.md)

| Type | Name |
| ---: | :--- |
|   | [**MatrixBatch**](classMatrixBatch.md#function-matrixbatch) (const std::size\_t batch\_size, const std::size\_t mat\_size) <br>_The constructor for_ [_**MatrixBatch**_](classMatrixBatch.md) _class._ |






## Detailed Description




**Template parameters:**


* `ExecSpace` Execution space,needed by Kokkos for allocations and parallelism. The simplest choice is to follow Kokkos, for that: specify Kokkos::DefaultExecutionSpace 
* `Solver` Refers to the solver type, default value is the Bicgstab which is more general. The use of a CG solver is also possible, in this case, please make sure that matrices structure fulfils CG requirements. 




    
## Public Functions Documentation




### function MatrixBatchCsr [1/2]

_The constructor for_ [_**MatrixBatchCsr**_](classMatrixBatchCsr.md) _class._
```C++
inline explicit MatrixBatchCsr::MatrixBatchCsr (
    const int batch_size,
    const int mat_size,
    const int nnz_per_system,
    std::optional< int > max_iter=std::nullopt,
    std::optional< double > res_tol=std::nullopt,
    std::optional< bool > logger=std::nullopt,
    std::optional< int > preconditionner_max_block_size=1u
) 
```





**Parameters:**


* `batch_size` Number of linear systems to solve. 
* `mat_size` Common matrix size for all the systems. 
* `nnz_per_system` Number of non-zero components per matrix. 
* `max_iter` maximal number of iterations for the solver, default 1000. 
* `res_tol` residual tolerance parameter, to ensure convergence. Be careful! the relative residual provided here, will be used as "implicit residual" in ginkgo solver. 
* `logger` boolean parameter for saving log information such residual and interactions count. 
* `preconditionner_max_block_size` An optional parameter used to define the maximum size of a block 




        

<hr>



### function MatrixBatchCsr [2/2]

_Constructor for_ [_**MatrixBatchCsr**_](classMatrixBatchCsr.md) _class._
```C++
inline explicit MatrixBatchCsr::MatrixBatchCsr (
    Kokkos::View< double **, Kokkos::LayoutRight, ExecSpace > batch_values,
    Kokkos::View< int *, Kokkos::LayoutRight, ExecSpace > cols_idx,
    Kokkos::View< int *, Kokkos::LayoutRight, ExecSpace > nnz_per_row,
    std::optional< int > max_iter=std::nullopt,
    std::optional< double > res_tol=std::nullopt,
    std::optional< bool > logger=std::nullopt,
    std::optional< int > preconditionner_max_block_size=1u
) 
```





**Parameters:**


* `batch_values` A 2D Kokkos view which stores the values of non-zero elements for the whole batch. 
* `cols_idx` A 1D Kokkos view which stores the column indices for each non-zero component.(only for one matrix) 
* `nnz_per_row` A 1D Kokkos view of length matrix\_size+1 which stores the count of the non-zeros along the lines of the matrix. It is defined as: nnz\_per\_row[0] = 0. nnz\_per\_row[matrix\_size] = total\_number\_of\_nonzero. To get the number of non-zero for a line i,one have to compute : n\_non\_zeros\_at\_line\_in = nnz\_per\_row[i+1]-nnz\_per\_row[i]. 
* `max_iter` maximal number of iterations for the solver, default 1000. 
* `res_tol` residual tolerance parameter, to ensure convergence. Be careful! The residual provided here, set as relative residual, will be used as "implicit residual" in ginkgo solver. Default value is set to 1e-15. 
* `logger` boolean parameter to save logger information. Default value false. 
* `preconditionner_max_block_size` An optional parameter used to define the maximum size of a block 




        

<hr>



### function get\_batch\_csr 

_A function to update information about values,indices and the number of non-zero per row for the whole batch. Data is managed by Kokkos Views stored on the host._ 
```C++
inline std::tuple< Kokkos::View< double **, Kokkos::LayoutRight, ExecSpace >, Kokkos::View< int *, Kokkos::LayoutRight, ExecSpace >, Kokkos::View< int *, Kokkos::LayoutRight, ExecSpace > > MatrixBatchCsr::get_batch_csr () 
```





**Returns:**

vals\_view A 2D Kokkos view which stores the values of non-zero elements for the whole batch. 




**Returns:**

col\_idx\_view A 1D Kokkos view which stores the column indices for each non-zero component.(only for one matrix) 




**Returns:**

nnz\_per\_row\_view A 1D Kokkos view which stores the count of non-zero per line, in an additive way. see nnz\_per\_row parameter in constructor. 





        

<hr>



### function norm 

_A function returning the norm of a matrix located at batch\_idx._ 
```C++
inline double MatrixBatchCsr::norm (
    int batch_idx
) const
```





**Parameters:**


* `batch_idx` The index of the matrix in the batch. 



**Returns:**

The value of the matrix infinite-norm. 





        

<hr>



### function setup\_solver 

_Perform a pre-process operation on the solver. Must be called after filling the matrix._ 
```C++
inline virtual void MatrixBatchCsr::setup_solver () 
```



It uses the batch of matrices to generate a batched Jacobi preconditioner. Other parameters like maximum number of iterations and tolerance are also used to instantiate a Ginkgo solver.


The stopping criterion is a reduction factor \|\|Ax-b\|\|/\|\|b\|\|&lt;tol with max\_iter maximum iterations. 


        
Implements [*MatrixBatch::setup\_solver*](classMatrixBatch.md#function-setup_solver)


<hr>



### function solve [1/2]

_Solve the batched linear problem Ax=b._ 
```C++
inline virtual void MatrixBatchCsr::solve (
    BatchedRHS const b
) const
```





**Parameters:**


* `b` A 2D Kokkos::View storing the batched right-hand sides of the problem and receiving the corresponding solutions. 




        
Implements [*MatrixBatch::solve*](classMatrixBatch.md#function-solve)


<hr>



### function solve [2/2]

_Solve the batched linear problem Ax=b._ 
```C++
inline void MatrixBatchCsr::solve (
    BatchedRHS const x,
    BatchedRHS const b
) const
```





**Parameters:**


* `x` A 2D Kokkos::View storing the batched initial guests (useful for iterative solver) of the problems, and receiving the corresponding solutions.
* `b` A 2D Kokkos::View storing the batched right-hand side of the problems. 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/matrix_tools/matrix_batch_csr.hpp`

