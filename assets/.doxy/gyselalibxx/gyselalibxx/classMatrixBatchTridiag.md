

# Class MatrixBatchTridiag

**template &lt;class ExecSpace&gt;**



[**ClassList**](annotated.md) **>** [**MatrixBatchTridiag**](classMatrixBatchTridiag.md)



_A structure for solving a set of independent tridiagonal systems using a direct method. The parallelism operates on the whole collection by dispatching to threads. Each problem is treated sequentially, by the tridiagonal matrix algorithm (TDMA). This solver is stable for tridiagonal matrices which satisfy one of the following conditions:_ [More...](#detailed-description)

* `#include <matrix_batch_tridiag.hpp>`



Inherits the following classes: [MatrixBatch](classMatrixBatch.md)
















## Public Types inherited from MatrixBatch

See [MatrixBatch](classMatrixBatch.md)

| Type | Name |
| ---: | :--- |
| typedef Kokkos::View&lt; double \*\*, Kokkos::LayoutRight, ExecSpace &gt; | [**BatchedRHS**](classMatrixBatch.md#typedef-batchedrhs)  <br>_The type of a Kokkos::View storing batched right-hand sides. Second dimension is batch dimension._  |






































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**MatrixBatchTridiag**](#function-matrixbatchtridiag) (const int batch\_size, const int mat\_size, DKokkosView2D const aa, DKokkosView2D const bb, DKokkosView2D const cc) <br>_Creates an instance of the_ [_**MatrixBatchTridiag**_](classMatrixBatchTridiag.md) _class. First dimension is the batch, second one refers to matrix entries indexed by line. The entries aa,bb,cc are 2D Kokkos views and have the same dimensions. LayoutRight: means that the "last" dimension is the contiguous one. aa(batch\_idx,0) and cc(batch\_idx,mat\_size-1) are not used for any values of batch\_idx._ |
|  bool | [**check\_stability**](#function-check_stability) () const<br>_Check if the matrices are in the stability area of the solver._  |
| virtual void | [**setup\_solver**](#function-setup_solver) () <br>_Perform a pre-process operation on the solver. Must be called after filling the matrix._  |
| virtual void | [**solve**](#function-solve) ([**BatchedRHS**](classMatrixBatch.md#typedef-batchedrhs) const b) const<br>_Solve the batched linear problem Ax=b._  |


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



* Diagonally Dominant.
* Symmetric positive-definite. Diagonally Dominant property is fully checked. Only symmetry property is checked, positivity-definiteness is not. 

**Template parameters:**


  * `ExecSpace` The execution space related to Kokkos. 






    
## Public Functions Documentation




### function MatrixBatchTridiag 

_Creates an instance of the_ [_**MatrixBatchTridiag**_](classMatrixBatchTridiag.md) _class. First dimension is the batch, second one refers to matrix entries indexed by line. The entries aa,bb,cc are 2D Kokkos views and have the same dimensions. LayoutRight: means that the "last" dimension is the contiguous one. aa(batch\_idx,0) and cc(batch\_idx,mat\_size-1) are not used for any values of batch\_idx._
```C++
inline explicit MatrixBatchTridiag::MatrixBatchTridiag (
    const int batch_size,
    const int mat_size,
    DKokkosView2D const aa,
    DKokkosView2D const bb,
    DKokkosView2D const cc
) 
```





**Parameters:**


* `batch_size` The size of the set of linear problems. 
* `mat_size` The common size of each individual matrix . 
* `aa` 2d Kokkos View which stores subdiagonal components for all matrices. 
* `bb` 2d Kokkos View which stores diagonal components for all matrices. 
* `cc` 2d Kokkos View which stores upper diagonal components for all matrices. 




        

<hr>



### function check\_stability 

_Check if the matrices are in the stability area of the solver._ 
```C++
inline bool MatrixBatchTridiag::check_stability () const
```



It checks if each matrix of the batch has one of the following structures:
* Diagonally Dominant.
* Symmetric. If assertion fails, there is at least one of the matrices which does not verify these conditions. Positivity-definiteness is not checked, it is needed to fully achieve TDMA algorithm requirement: symmetric positive definite.






**Returns:**

A boolean which indicates if the stability condition may be verified. 





        

<hr>



### function setup\_solver 

_Perform a pre-process operation on the solver. Must be called after filling the matrix._ 
```C++
inline virtual void MatrixBatchTridiag::setup_solver () 
```



It calls check\_stability function to verify if the matrices data is in range of validity of the solver.


The stopping criterion is a reduction factor \|\|Ax-b\|\|/\|\|b\|\|&lt;tol with max\_iter maximum iterations. 


        
Implements [*MatrixBatch::setup\_solver*](classMatrixBatch.md#function-setup_solver)


<hr>



### function solve 

_Solve the batched linear problem Ax=b._ 
```C++
inline virtual void MatrixBatchTridiag::solve (
    BatchedRHS const b
) const
```





**Parameters:**


* `b` A 2D Kokkos::View storing the batched right-hand sides of the problem and receiving the corresponding solutions. 




        
Implements [*MatrixBatch::solve*](classMatrixBatch.md#function-solve)


<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/matrix_tools/matrix_batch_tridiag.hpp`

