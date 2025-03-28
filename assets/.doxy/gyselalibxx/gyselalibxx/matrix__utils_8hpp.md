

# File matrix\_utils.hpp



[**FileList**](files.md) **>** [**matrix\_tools**](dir_8cedd1260cc2f2819c8df2fc66ad98b5.md) **>** [**matrix\_utils.hpp**](matrix__utils_8hpp.md)

[Go to the source code of this file](matrix__utils_8hpp_source.md)



* `#include <ginkgo/ginkgo.hpp>`
* `#include <Kokkos_Core.hpp>`





































## Public Functions

| Type | Name |
| ---: | :--- |
|  void | [**check\_conv**](#function-check_conv) (int const batch\_size, double const tol, std::shared\_ptr&lt; const gko::Executor &gt; gko\_exec, std::shared\_ptr&lt; const gko::batch::log::BatchConvergence&lt; double &gt; &gt; logger) <br>_A function for checking convergence. It loops over the batch and checks the if residual is lower or equal to the prescribed tolerance._  |
|  unsigned int | [**default\_preconditionner\_max\_block\_size**](#function-default_preconditionner_max_block_size) () noexcept<br> |
|  void | [**save\_logger**](#function-save_logger) (std::fstream & log\_file, int const batch\_index, std::unique\_ptr&lt; sparse\_type &gt; matrix, Kokkos::View&lt; double \*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace &gt; const x\_view, Kokkos::View&lt; double \*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace &gt; const b\_view, std::shared\_ptr&lt; const gko::log::Convergence&lt; double &gt; &gt; logger, double const tol) <br>_A function to save convergence data using the logger._  |
|  void | [**save\_logger**](#function-save_logger) (std::fstream & log\_file, std::shared\_ptr&lt; batch\_sparse\_type &gt; batch\_matrix, Kokkos::View&lt; double \*\*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace &gt; const x\_view, Kokkos::View&lt; double \*\*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace &gt; const b\_view, std::shared\_ptr&lt; const gko::batch::log::BatchConvergence&lt; double &gt; &gt; logger, double const tol) <br>_A function to save convergence data using the logger._  |
|  auto | [**to\_gko\_multivector**](#function-to_gko_multivector) (std::shared\_ptr&lt; const gko::Executor &gt; const & gko\_exec, KokkosViewType const & view) <br>_A function to convert a 2D Kokkos view into a ginkgo multivector structure._  |
|  void | [**write\_log**](#function-write_log) (std::fstream & log\_file, int const batch\_index, int const num\_iterations, double const implicit\_res\_norm, double const true\_res\_norm, double const b\_norm, double const tol) <br>_A helper to write the log corresponding to a single batch._  |




























## Public Functions Documentation




### function check\_conv 

_A function for checking convergence. It loops over the batch and checks the if residual is lower or equal to the prescribed tolerance._ 
```C++
inline void check_conv (
    int const batch_size,
    double const tol,
    std::shared_ptr< const gko::Executor > gko_exec,
    std::shared_ptr< const gko::batch::log::BatchConvergence< double > > logger
) 
```





**Parameters:**


* `batch_size` the size of the batch , ie number of linears problems. 
* `tol` The tolerancy on residual norm above which a non-convergency is reported. 
* `gko_exec` Ginkgo executor, refers to the execution space. 
* `logger` Ginkgo convergence object which stores iterations number and residual for the whole batch. 




        

<hr>



### function default\_preconditionner\_max\_block\_size 

```C++
template<class ExecSpace>
unsigned int default_preconditionner_max_block_size () noexcept
```




<hr>



### function save\_logger 

_A function to save convergence data using the logger._ 
```C++
template<class sparse_type>
void save_logger (
    std::fstream & log_file,
    int const batch_index,
    std::unique_ptr< sparse_type > matrix,
    Kokkos::View< double *, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace > const x_view,
    Kokkos::View< double *, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace > const b_view,
    std::shared_ptr< const gko::log::Convergence< double > > logger,
    double const tol
) 
```





**Parameters:**


* `log_file` The file in which the logs will be saved. 
* `batch_index` The index of the relevant matrix in the batch of matrices. 
* `matrix` Batch of matrices. 
* `x_view` 2d Kokkos view containing the batch of computed solutions. 
* `b_view` 2d Kokkos view containing the batch of rhs. 
* `logger` Ginkgo logger which stores residual and numbers of iterations for the whole batch. 
* `tol` tolerance 




        

<hr>



### function save\_logger 

_A function to save convergence data using the logger._ 
```C++
template<class batch_sparse_type>
void save_logger (
    std::fstream & log_file,
    std::shared_ptr< batch_sparse_type > batch_matrix,
    Kokkos::View< double **, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace > const x_view,
    Kokkos::View< double **, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace > const b_view,
    std::shared_ptr< const gko::batch::log::BatchConvergence< double > > logger,
    double const tol
) 
```





**Parameters:**


* `log_file` The file in which the logs will be saved. 
* `batch_matrix` Batch of matrices. 
* `x_view` 2d Kokkos view containing the batch of computed solutions. 
* `b_view` 2d Kokkos view containing the batch of rhs. 
* `logger` Ginkgo logger which stores residual and numbers of iterations for the whole batch. 
* `tol` tolerance 




        

<hr>



### function to\_gko\_multivector 

_A function to convert a 2D Kokkos view into a ginkgo multivector structure._ 
```C++
template<class KokkosViewType>
auto to_gko_multivector (
    std::shared_ptr< const gko::Executor > const & gko_exec,
    KokkosViewType const & view
) 
```





**Parameters:**


* `gko_exec` A Ginkgo executor that has access to the Kokkos::View memory space 
* `view` A 2-D Kokkos::View with unit stride in the second dimension 



**Returns:**

A Ginkgo Multivector view over the Kokkos::View data 





        

<hr>



### function write\_log 

_A helper to write the log corresponding to a single batch._ 
```C++
inline void write_log (
    std::fstream & log_file,
    int const batch_index,
    int const num_iterations,
    double const implicit_res_norm,
    double const true_res_norm,
    double const b_norm,
    double const tol
) 
```





**Parameters:**


* `log_file` The stream of the log file. 
* `batch_index` The index of the batch. 
* `num_iterations` The number of iterations the iterative solver performed for this batch. 
* `implicit_res_norm` The implicit residual norm at the end of the solver call (evaluated by the Ginkgo solver). 
* `true_res_norm` The true residual norm at the end of the solver call (re-computed by hand). 
* `b_norm` The norm of the right-hand side. 
* `tol` The tolerancy on residual norm above which a non-convergency is reported. 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/matrix_tools/matrix_utils.hpp`

