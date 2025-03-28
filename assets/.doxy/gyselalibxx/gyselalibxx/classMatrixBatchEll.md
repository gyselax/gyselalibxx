

# Class MatrixBatchEll

**template &lt;class ExecSpace&gt;**



[**ClassList**](annotated.md) **>** [**MatrixBatchEll**](classMatrixBatchEll.md)



[_**Matrix**_](classMatrix.md) _class which is able to manage and solve a batch of sparse linear systems. Executes on either CPU or GPU. It takes advantage of the sparse structure, and the only batched solver available in Ginkgo : Stabilised Bicg. The sparsity pattern is assumed to be the same for all matrices. ie the non-zero components are located at the same places for all matrices. This class uses the ELL storage format which needs two 1D arrays, one stores values the other column indices. The class returns these arrays (as Kokkos views) with the get\_batch\_idx\_and\_vals function, it is then possible to fill them outside the class. Tolerance and maximal number of iterations, which are parameters for the iterative solver, are set in the constructor. It is possible to get convergence information by activating the logger at constructor call._[More...](#detailed-description)

* `#include <matrix_batch_ell.hpp>`



Inherits the following classes: [MatrixBatch](classMatrixBatch.md)
















## Public Types inherited from MatrixBatch

See [MatrixBatch](classMatrixBatch.md)

| Type | Name |
| ---: | :--- |
| typedef Kokkos::View&lt; double \*\*, Kokkos::LayoutRight, ExecSpace &gt; | [**BatchedRHS**](classMatrixBatch.md#typedef-batchedrhs)  <br>_The type of a Kokkos::View storing batched right-hand sides. Second dimension is batch dimension._  |






































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**MatrixBatchEll**](#function-matrixbatchell-12) (const int batch\_size, const int mat\_size, const int non\_zeros\_per\_row, std::optional&lt; int &gt; max\_iter=std::nullopt, std::optional&lt; double &gt; res\_tol=std::nullopt, std::optional&lt; bool &gt; logger=std::nullopt) <br>_The constructor for_ [_**MatrixBatchEll**_](classMatrixBatchEll.md) _class._ |
|   | [**MatrixBatchEll**](#function-matrixbatchell-22) (Kokkos::View&lt; int \*\*, Kokkos::LayoutLeft, ExecSpace &gt; cols\_idx, Kokkos::View&lt; double \*\*\*, Kokkos::LayoutStride, ExecSpace &gt; batch\_values, std::optional&lt; int &gt; max\_iter=std::nullopt, std::optional&lt; double &gt; res\_tol=std::nullopt, std::optional&lt; bool &gt; logger=std::nullopt) <br>_Constructor for_ [_**MatrixBatchEll**_](classMatrixBatchEll.md) _class._ |
|  std::pair&lt; Kokkos::View&lt; int \*\*, Kokkos::LayoutLeft, ExecSpace &gt;, Kokkos::View&lt; double \*\*\*, Kokkos::LayoutStride, ExecSpace &gt; &gt; | [**get\_batch\_ell**](#function-get_batch_ell) () <br>_A function to get information about values and indices for the whole batch. Data is managed by two Kokkos Views stored on the host._  |
|  double | [**get\_ell\_element**](#function-get_ell_element) (int batch\_idx, int line\_idx, int non\_zero\_col\_idx) const<br>_A getter function for a value located at a specified place._  |
|  double | [**norm**](#function-norm) (int batch\_idx) const<br>_A function returns the norm of a matrix located at batch\_idx._  |
|  void | [**set\_ell\_element**](#function-set_ell_element) (int batch\_idx, int line\_idx, int non\_zero\_col\_idx, double aij) <br>_A setter function to modify a value located at a specified place._  |
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




**Template parameters:**


* `ExecSpace` Execution space,needed by Kokkos for allocations and parallelism. The simplest choice is to follow Kokkos, for that: specify Kokkos::DefaultExecutionSpace 




    
## Public Functions Documentation




### function MatrixBatchEll [1/2]

_The constructor for_ [_**MatrixBatchEll**_](classMatrixBatchEll.md) _class._
```C++
inline explicit MatrixBatchEll::MatrixBatchEll (
    const int batch_size,
    const int mat_size,
    const int non_zeros_per_row,
    std::optional< int > max_iter=std::nullopt,
    std::optional< double > res_tol=std::nullopt,
    std::optional< bool > logger=std::nullopt
) 
```





**Parameters:**


* `batch_size` Number of linear systems to solve. 
* `mat_size` Common matrix size for all the systems. 
* `non_zeros_per_row` number of non zero components per line. 
* `max_iter` maximal number of iterations for the solver 
* `res_tol` residual tolerance parameter, to ensure convergence. Be careful! the relative residual provided here, will be used as "implicit residual" in ginkgo solver. 
* `logger` boolean parameter for saving log information such residual and interactions count. 




        

<hr>



### function MatrixBatchEll [2/2]

_Constructor for_ [_**MatrixBatchEll**_](classMatrixBatchEll.md) _class._
```C++
inline explicit MatrixBatchEll::MatrixBatchEll (
    Kokkos::View< int **, Kokkos::LayoutLeft, ExecSpace > cols_idx,
    Kokkos::View< double ***, Kokkos::LayoutStride, ExecSpace > batch_values,
    std::optional< int > max_iter=std::nullopt,
    std::optional< double > res_tol=std::nullopt,
    std::optional< bool > logger=std::nullopt
) 
```





**Parameters:**


* `cols_idx` A Kokkos view which stores the column indices of non-zero components. 
* `batch_values` A Kokkos view which stores the values of non-zero elements. 
* `max_iter` maximal number of iterations for the solver, default 500. 
* `res_tol` residual tolerance parameter, to ensure convergence. Be careful! The residual provided here, set as relative residual, will be used as "implicit residual" in ginkgo solver. Default value is set to 1e-15. 
* `logger` boolean parameter to save logger information. Default value false. 




        

<hr>



### function get\_batch\_ell 

_A function to get information about values and indices for the whole batch. Data is managed by two Kokkos Views stored on the host._ 
```C++
inline std::pair< Kokkos::View< int **, Kokkos::LayoutLeft, ExecSpace >, Kokkos::View< double ***, Kokkos::LayoutStride, ExecSpace > > MatrixBatchEll::get_batch_ell () 
```





**Returns:**

idx\_view Column indices for the non-zero values. 




**Returns:**

vals\_view The non-zero values. 





        

<hr>



### function get\_ell\_element 

_A getter function for a value located at a specified place._ 
```C++
inline double MatrixBatchEll::get_ell_element (
    int batch_idx,
    int line_idx,
    int non_zero_col_idx
) const
```





**Parameters:**


* `batch_idx` Index in the batch. 
* `line_idx` Line index inside the matrix. 
* `non_zero_col_idx` Non-zero index element in the line. 



**Returns:**

value of the component. 





        

<hr>



### function norm 

_A function returns the norm of a matrix located at batch\_idx._ 
```C++
inline double MatrixBatchEll::norm (
    int batch_idx
) const
```





**Parameters:**


* `batch_idx` integer, index of the matrix in the batch. 



**Returns:**

value of the matrix infinite-norm. 





        

<hr>



### function set\_ell\_element 

_A setter function to modify a value located at a specified place._ 
```C++
inline void MatrixBatchEll::set_ell_element (
    int batch_idx,
    int line_idx,
    int non_zero_col_idx,
    double aij
) 
```





**Parameters:**


* `batch_idx` Index in the batch. 
* `line_idx` Line index inside the matrix. 
* `non_zero_col_idx` Non-zero index element in the line. 
* `aij` New value. 




        

<hr>



### function setup\_solver 

_Perform a pre-process operation on the solver. Must be called after filling the matrix._ 
```C++
inline virtual void MatrixBatchEll::setup_solver () 
```



It uses parameters like maximum number of iterations and tolerance are used to instantiate a Ginkgo solver.


The stopping criterion is a reduction factor \|\|Ax-b\|\|/\|\|b\|\|&lt;tol with max\_iter maximum iterations. 


        
Implements [*MatrixBatch::setup\_solver*](classMatrixBatch.md#function-setup_solver)


<hr>



### function solve 

_Solve the batched linear problem Ax=b._ 
```C++
inline virtual void MatrixBatchEll::solve (
    BatchedRHS const b
) const
```





**Parameters:**


* `b` A 2D Kokkos::View storing the batched right-hand sides of the problem and receiving the corresponding solutions. 




        
Implements [*MatrixBatch::solve*](classMatrixBatch.md#function-solve)


<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/matrix_tools/matrix_batch_ell.hpp`

