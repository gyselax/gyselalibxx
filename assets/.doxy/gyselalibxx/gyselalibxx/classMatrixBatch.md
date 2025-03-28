

# Class MatrixBatch

**template &lt;typename ExecSpace&gt;**



[**ClassList**](annotated.md) **>** [**MatrixBatch**](classMatrixBatch.md)



[_**MatrixBatch**_](classMatrixBatch.md) _superclass for managing a collection of linear systems. The main assumption is that all matrices have the same size. It is also assumed that each matrix is used to solve one equation._[More...](#detailed-description)

* `#include <matrix_batch.hpp>`





Inherited by the following classes: [MatrixBatchCsr](classMatrixBatchCsr.md),  [MatrixBatchCsr](classMatrixBatchCsr.md),  [MatrixBatchEll](classMatrixBatchEll.md),  [MatrixBatchTridiag](classMatrixBatchTridiag.md)












## Public Types

| Type | Name |
| ---: | :--- |
| typedef Kokkos::View&lt; double \*\*, Kokkos::LayoutRight, ExecSpace &gt; | [**BatchedRHS**](#typedef-batchedrhs)  <br>_The type of a Kokkos::View storing batched right-hand sides. Second dimension is batch dimension._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
|  std::size\_t | [**batch\_size**](#function-batch_size) () const<br>_Get the batch size of the linear problem._  |
| virtual void | [**setup\_solver**](#function-setup_solver) () = 0<br>_Perform a pre-process operation on the solver. Must be called after filling the matrix._  |
|  std::size\_t | [**size**](#function-size) () const<br>_Get the size of the square matrix corresponding to a single batch in one of its dimensions._  |
| virtual void | [**solve**](#function-solve) ([**BatchedRHS**](classMatrixBatch.md#typedef-batchedrhs) b) const = 0<br>_Solve the multiple right-hand sides linear problem Ax=b inplace._  |
| virtual  | [**~MatrixBatch**](#function-matrixbatch) () = default<br>_Destruct._  |
























## Protected Functions

| Type | Name |
| ---: | :--- |
|   | [**MatrixBatch**](#function-matrixbatch) (const std::size\_t batch\_size, const std::size\_t mat\_size) <br>_The constructor for_ [_**MatrixBatch**_](classMatrixBatch.md) _class._ |




## Detailed Description


Classes inheriting from this class must manage other aspects: Sparsity: kind of storage (Dense, Ell, Csr, etc.) Kind of solver (direct, iterative) Preconditioners and factorisations




**Template parameters:**


* `ExecSpace` Execution space,needed by Kokkos for allocations and parallelism. The simplest choice is to follow Kokkos, for that: specify Kokkos::DefaultExecutionSpace 




    
## Public Types Documentation




### typedef BatchedRHS 

_The type of a Kokkos::View storing batched right-hand sides. Second dimension is batch dimension._ 
```C++
using MatrixBatch< ExecSpace >::BatchedRHS =  Kokkos::View<double**, Kokkos::LayoutRight, ExecSpace>;
```




<hr>
## Public Functions Documentation




### function batch\_size 

_Get the batch size of the linear problem._ 
```C++
inline std::size_t MatrixBatch::batch_size () const
```





**Returns:**

The batch size of the linear problem. 





        

<hr>



### function setup\_solver 

_Perform a pre-process operation on the solver. Must be called after filling the matrix._ 
```C++
virtual void MatrixBatch::setup_solver () = 0
```




<hr>



### function size 

_Get the size of the square matrix corresponding to a single batch in one of its dimensions._ 
```C++
inline std::size_t MatrixBatch::size () const
```





**Returns:**

The size of the matrix in one of its dimensions. 





        

<hr>



### function solve 

_Solve the multiple right-hand sides linear problem Ax=b inplace._ 
```C++
virtual void MatrixBatch::solve (
    BatchedRHS b
) const = 0
```





**Parameters:**


* `b` A 2D Kokkos::View storing the batched right-hand sides of the problem and receiving the corresponding solution. 




        

<hr>



### function ~MatrixBatch 

_Destruct._ 
```C++
virtual MatrixBatch::~MatrixBatch () = default
```




<hr>
## Protected Functions Documentation




### function MatrixBatch 

_The constructor for_ [_**MatrixBatch**_](classMatrixBatch.md) _class._
```C++
inline explicit MatrixBatch::MatrixBatch (
    const std::size_t batch_size,
    const std::size_t mat_size
) 
```





**Parameters:**


* `batch_size` Number of linear systems to solve. 
* `mat_size` Common matrix size for all the systems. 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/matrix_tools/matrix_batch.hpp`

