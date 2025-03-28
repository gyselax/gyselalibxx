

# Class FFTPoissonSolver&lt; IdxRange&lt; GridPDEDim1D... &gt;, IdxRangeFull, ExecSpace, LayoutSpace &gt;

**template &lt;class... GridPDEDim1D, class IdxRangeFull, class ExecSpace, class LayoutSpace&gt;**



[**ClassList**](annotated.md) **>** [**FFTPoissonSolver&lt; IdxRange&lt; GridPDEDim1D... &gt;, IdxRangeFull, ExecSpace, LayoutSpace &gt;**](classFFTPoissonSolver_3_01IdxRange_3_01GridPDEDim1D_8_8_8_01_4_00_01IdxRangeFull_00_01ExecSpace_00_01LayoutSpace_01_4.md)



_A class to solve the following equation:_  _using a Fourier transform._[More...](#detailed-description)

* `#include <fft_poisson_solver.hpp>`



Inherits the following classes: [IPoissonSolver](classIPoissonSolver.md)












## Classes

| Type | Name |
| ---: | :--- |
| struct | [**GridFourier**](structFFTPoissonSolver_3_01IdxRange_3_01GridPDEDim1D_8_8_8_01_4_00_01IdxRangeFull_00_01ExecSpace2aeecfe91d464f5738599cc105fb6087.md) &lt;class Dim&gt;<br> |


## Public Types

| Type | Name |
| ---: | :--- |
| typedef typename base\_type::batch\_idx\_range\_type | [**batch\_idx\_range\_type**](#typedef-batch_idx_range_type)  <br>_The index range type describing the batch dimensions._  |
| typedef typename base\_type::batch\_index\_type | [**batch\_index\_type**](#typedef-batch_index_type)  <br>_The index type for indexing a batch dimension._  |
| typedef typename base\_type::const\_field\_type | [**const\_field\_type**](#typedef-const_field_type)  <br>_The const Field type of the arguments to operator()._  |
| typedef typename base\_type::field\_type | [**field\_type**](#typedef-field_type)  <br>_The Field type of the arguments to operator()._  |
| typedef FieldMem&lt; Kokkos::complex&lt; double &gt;, [**fourier\_idx\_range\_type**](classFFTPoissonSolver_3_01IdxRange_3_01GridPDEDim1D_8_8_8_01_4_00_01IdxRangeFull_00_01ExecSpace_00_01LayoutSpace_01_4.md#typedef-fourier_idx_range_type), [**memory\_space**](classFFTPoissonSolver_3_01IdxRange_3_01GridPDEDim1D_8_8_8_01_4_00_01IdxRangeFull_00_01ExecSpace_00_01LayoutSpace_01_4.md#typedef-memory_space) &gt; | [**fourier\_field\_mem\_type**](#typedef-fourier_field_mem_type)  <br>_The type of a Field storing the Fourier transform of a function._  |
| typedef typename fourier\_field\_mem\_type::span\_type | [**fourier\_field\_type**](#typedef-fourier_field_type)  <br>_The type of a Field storing the Fourier transform of a function._  |
| typedef IdxRange&lt; GridFourier&lt; typename GridPDEDim1D::continuous\_dimension\_type &gt;... &gt; | [**fourier\_idx\_range\_type**](#typedef-fourier_idx_range_type)  <br>_The type of the Fourier space index range._  |
| typedef typename fourier\_idx\_range\_type::discrete\_element\_type | [**fourier\_index\_type**](#typedef-fourier_index_type)  <br>_The type of an index of the Fourier space index range._  |
| typedef typename base\_type::laplacian\_idx\_range\_type | [**laplacian\_idx\_range\_type**](#typedef-laplacian_idx_range_type)  <br>_The type of the index range on which the equation is defined._  |
| typedef typename base\_type::layout\_space | [**layout\_space**](#typedef-layout_space)  <br>_The layout space of the Fields passed to operator()._  |
| typedef typename base\_type::memory\_space | [**memory\_space**](#typedef-memory_space)  <br>_The space (CPU/GPU) where the Fields passed to operator() are saved._  |
| typedef typename base\_type::vector\_field\_type | [**vector\_field\_type**](#typedef-vector_field_type)  <br>_The type of the derivative of_  _._ |








































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**FFTPoissonSolver**](#function-fftpoissonsolver) ([**laplacian\_idx\_range\_type**](classFFTPoissonSolver_3_01IdxRange_3_01GridPDEDim1D_8_8_8_01_4_00_01IdxRangeFull_00_01ExecSpace_00_01LayoutSpace_01_4.md#typedef-laplacian_idx_range_type) laplacian\_idx\_range) <br>_A constructor for the FFT Poisson solver. This constructor calls ddc::init\_discrete\_space so it should only be called once per simulation._  |
|  void | [**negative\_differentiate\_equation**](#function-negative_differentiate_equation) ([**fourier\_field\_type**](classFFTPoissonSolver_3_01IdxRange_3_01GridPDEDim1D_8_8_8_01_4_00_01IdxRangeFull_00_01ExecSpace_00_01LayoutSpace_01_4.md#typedef-fourier_field_type) derivative, [**fourier\_field\_type**](classFFTPoissonSolver_3_01IdxRange_3_01GridPDEDim1D_8_8_8_01_4_00_01IdxRangeFull_00_01ExecSpace_00_01LayoutSpace_01_4.md#typedef-fourier_field_type) values) const<br>_Differentiate and multiply by -1 an expression in Fourier space by multiplying by -i \* k This function should be private. It is not due to the inclusion of a KOKKOS\_LAMBDA._  |
| virtual [**field\_type**](classFFTPoissonSolver_3_01IdxRange_3_01GridPDEDim1D_8_8_8_01_4_00_01IdxRangeFull_00_01ExecSpace_00_01LayoutSpace_01_4.md#typedef-field_type) | [**operator()**](#function-operator) ([**field\_type**](classFFTPoissonSolver_3_01IdxRange_3_01GridPDEDim1D_8_8_8_01_4_00_01IdxRangeFull_00_01ExecSpace_00_01LayoutSpace_01_4.md#typedef-field_type) phi, [**field\_type**](classFFTPoissonSolver_3_01IdxRange_3_01GridPDEDim1D_8_8_8_01_4_00_01IdxRangeFull_00_01ExecSpace_00_01LayoutSpace_01_4.md#typedef-field_type) rho) const<br>_An operator which calculates the solution_  _to Poisson's equation:_ _._ |
| virtual [**field\_type**](classFFTPoissonSolver_3_01IdxRange_3_01GridPDEDim1D_8_8_8_01_4_00_01IdxRangeFull_00_01ExecSpace_00_01LayoutSpace_01_4.md#typedef-field_type) | [**operator()**](#function-operator_1) ([**field\_type**](classFFTPoissonSolver_3_01IdxRange_3_01GridPDEDim1D_8_8_8_01_4_00_01IdxRangeFull_00_01ExecSpace_00_01LayoutSpace_01_4.md#typedef-field_type) phi, [**vector\_field\_type**](classFFTPoissonSolver_3_01IdxRange_3_01GridPDEDim1D_8_8_8_01_4_00_01IdxRangeFull_00_01ExecSpace_00_01LayoutSpace_01_4.md#typedef-vector_field_type) E, [**field\_type**](classFFTPoissonSolver_3_01IdxRange_3_01GridPDEDim1D_8_8_8_01_4_00_01IdxRangeFull_00_01ExecSpace_00_01LayoutSpace_01_4.md#typedef-field_type) rho) const<br>_An operator which calculates the solution_  _to Poisson's equation and its derivative:_ __ _._ |
|  void | [**solve\_poisson\_equation**](#function-solve_poisson_equation) ([**fourier\_field\_type**](classFFTPoissonSolver_3_01IdxRange_3_01GridPDEDim1D_8_8_8_01_4_00_01IdxRangeFull_00_01ExecSpace_00_01LayoutSpace_01_4.md#typedef-fourier_field_type) intermediate\_chunk, DField&lt; [**laplacian\_idx\_range\_type**](classFFTPoissonSolver_3_01IdxRange_3_01GridPDEDim1D_8_8_8_01_4_00_01IdxRangeFull_00_01ExecSpace_00_01LayoutSpace_01_4.md#typedef-laplacian_idx_range_type), [**memory\_space**](classFFTPoissonSolver_3_01IdxRange_3_01GridPDEDim1D_8_8_8_01_4_00_01IdxRangeFull_00_01ExecSpace_00_01LayoutSpace_01_4.md#typedef-memory_space), Layout &gt; rho) const<br>_A function to solve the Poisson equation in Fourier space This function should be private. It is not due to the inclusion of a KOKKOS\_LAMBDA._  |
























































## Detailed Description


The implementation of this class can be found at [**FFTPoissonSolver**](classFFTPoissonSolver.md)&lt; IdxRange&lt;GridPDEDim1D...&gt;, IdxRangeFull, ExecSpace, LayoutSpace &gt;.  

**Template parameters:**


* `IdxRangeLaplacian` The index range on which the equation is defined. 
* `IdxRangeFull` The index range on which the operator() acts. This is equal to the IdxRangeLaplacian plus any batched dimensions. 
* `ExecSpace` The space (CPU/GPU) where the calculations will take place. 
* `LayoutSpace` The layout space of the Fields passed to operator(). 




    
## Public Types Documentation




### typedef batch\_idx\_range\_type 

_The index range type describing the batch dimensions._ 
```C++
using FFTPoissonSolver< IdxRange< GridPDEDim1D... >, IdxRangeFull, ExecSpace, LayoutSpace >::batch_idx_range_type =  typename base_type::batch_idx_range_type;
```




<hr>



### typedef batch\_index\_type 

_The index type for indexing a batch dimension._ 
```C++
using FFTPoissonSolver< IdxRange< GridPDEDim1D... >, IdxRangeFull, ExecSpace, LayoutSpace >::batch_index_type =  typename base_type::batch_index_type;
```




<hr>



### typedef const\_field\_type 

_The const Field type of the arguments to operator()._ 
```C++
using FFTPoissonSolver< IdxRange< GridPDEDim1D... >, IdxRangeFull, ExecSpace, LayoutSpace >::const_field_type =  typename base_type::const_field_type;
```




<hr>



### typedef field\_type 

_The Field type of the arguments to operator()._ 
```C++
using FFTPoissonSolver< IdxRange< GridPDEDim1D... >, IdxRangeFull, ExecSpace, LayoutSpace >::field_type =  typename base_type::field_type;
```




<hr>



### typedef fourier\_field\_mem\_type 

_The type of a Field storing the Fourier transform of a function._ 
```C++
using FFTPoissonSolver< IdxRange< GridPDEDim1D... >, IdxRangeFull, ExecSpace, LayoutSpace >::fourier_field_mem_type =  FieldMem<Kokkos::complex<double>, fourier_idx_range_type, memory_space>;
```




<hr>



### typedef fourier\_field\_type 

_The type of a Field storing the Fourier transform of a function._ 
```C++
using FFTPoissonSolver< IdxRange< GridPDEDim1D... >, IdxRangeFull, ExecSpace, LayoutSpace >::fourier_field_type =  typename fourier_field_mem_type::span_type;
```




<hr>



### typedef fourier\_idx\_range\_type 

_The type of the Fourier space index range._ 
```C++
using FFTPoissonSolver< IdxRange< GridPDEDim1D... >, IdxRangeFull, ExecSpace, LayoutSpace >::fourier_idx_range_type =  IdxRange<GridFourier<typename GridPDEDim1D::continuous_dimension_type>...>;
```




<hr>



### typedef fourier\_index\_type 

_The type of an index of the Fourier space index range._ 
```C++
using FFTPoissonSolver< IdxRange< GridPDEDim1D... >, IdxRangeFull, ExecSpace, LayoutSpace >::fourier_index_type =  typename fourier_idx_range_type::discrete_element_type;
```




<hr>



### typedef laplacian\_idx\_range\_type 

_The type of the index range on which the equation is defined._ 
```C++
using FFTPoissonSolver< IdxRange< GridPDEDim1D... >, IdxRangeFull, ExecSpace, LayoutSpace >::laplacian_idx_range_type =  typename base_type::laplacian_idx_range_type;
```




<hr>



### typedef layout\_space 

_The layout space of the Fields passed to operator()._ 
```C++
using FFTPoissonSolver< IdxRange< GridPDEDim1D... >, IdxRangeFull, ExecSpace, LayoutSpace >::layout_space =  typename base_type::layout_space;
```




<hr>



### typedef memory\_space 

_The space (CPU/GPU) where the Fields passed to operator() are saved._ 
```C++
using FFTPoissonSolver< IdxRange< GridPDEDim1D... >, IdxRangeFull, ExecSpace, LayoutSpace >::memory_space =  typename base_type::memory_space;
```




<hr>



### typedef vector\_field\_type 

_The type of the derivative of_  _._
```C++
using FFTPoissonSolver< IdxRange< GridPDEDim1D... >, IdxRangeFull, ExecSpace, LayoutSpace >::vector_field_type =  typename base_type::vector_field_type;
```




<hr>
## Public Functions Documentation




### function FFTPoissonSolver 

_A constructor for the FFT Poisson solver. This constructor calls ddc::init\_discrete\_space so it should only be called once per simulation._ 
```C++
inline explicit FFTPoissonSolver< IdxRange< GridPDEDim1D... >, IdxRangeFull, ExecSpace, LayoutSpace >::FFTPoissonSolver (
    laplacian_idx_range_type laplacian_idx_range
) 
```





**Parameters:**


* `laplacian_idx_range` The index range on which the equation should be solved. 




        

<hr>



### function negative\_differentiate\_equation 

_Differentiate and multiply by -1 an expression in Fourier space by multiplying by -i \* k This function should be private. It is not due to the inclusion of a KOKKOS\_LAMBDA._ 
```C++
template<class Dim>
inline void FFTPoissonSolver< IdxRange< GridPDEDim1D... >, IdxRangeFull, ExecSpace, LayoutSpace >::negative_differentiate_equation (
    fourier_field_type derivative,
    fourier_field_type values
) const
```





**Parameters:**


* `derivative` The Field where the derivative will be saved. 
* `values` The Field containing the values of the function in Fourier space.



**Template parameters:**


* `Dim` The dimension along which the expression is differentiated. 




        

<hr>



### function operator() 

_An operator which calculates the solution_  _to Poisson's equation:_ _._
```C++
inline virtual field_type FFTPoissonSolver< IdxRange< GridPDEDim1D... >, IdxRangeFull, ExecSpace, LayoutSpace >::operator() (
    field_type phi,
    field_type rho
) const
```





**Parameters:**


* `phi` The solution to Poisson's equation. 
* `rho` The right-hand side of Poisson's equation.



**Returns:**

A reference to the solution to Poisson's equation. 





        

<hr>



### function operator() 

_An operator which calculates the solution_  _to Poisson's equation and its derivative:_ __ _._
```C++
inline virtual field_type FFTPoissonSolver< IdxRange< GridPDEDim1D... >, IdxRangeFull, ExecSpace, LayoutSpace >::operator() (
    field_type phi,
    vector_field_type E,
    field_type rho
) const
```





**Parameters:**


* `phi` The solution to Poisson's equation. 
* `E` The derivative of the solution to Poisson's equation. 
* `rho` The right-hand side of Poisson's equation.



**Returns:**

A reference to the solution to Poisson's equation. 





        

<hr>



### function solve\_poisson\_equation 

_A function to solve the Poisson equation in Fourier space This function should be private. It is not due to the inclusion of a KOKKOS\_LAMBDA._ 
```C++
template<class Layout>
inline void FFTPoissonSolver< IdxRange< GridPDEDim1D... >, IdxRangeFull, ExecSpace, LayoutSpace >::solve_poisson_equation (
    fourier_field_type intermediate_chunk,
    DField< laplacian_idx_range_type , memory_space , Layout > rho
) const
```





**Parameters:**


* `intermediate_chunk` The solution to the Poisson equation in Fourier space. 
* `rho` The right-hand side of the Poisson equation in Fourier space. 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/pde_solvers/fft_poisson_solver.hpp`

