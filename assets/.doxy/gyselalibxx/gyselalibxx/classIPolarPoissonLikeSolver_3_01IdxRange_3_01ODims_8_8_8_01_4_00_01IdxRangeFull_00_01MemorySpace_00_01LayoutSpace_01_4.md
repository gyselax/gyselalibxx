

# Class IPolarPoissonLikeSolver&lt; IdxRange&lt; ODims... &gt;, IdxRangeFull, MemorySpace, LayoutSpace &gt;

**template &lt;class... ODims, class IdxRangeFull, class MemorySpace, class LayoutSpace&gt;**



[**ClassList**](annotated.md) **>** [**IPolarPoissonLikeSolver&lt; IdxRange&lt; ODims... &gt;, IdxRangeFull, MemorySpace, LayoutSpace &gt;**](classIPolarPoissonLikeSolver_3_01IdxRange_3_01ODims_8_8_8_01_4_00_01IdxRangeFull_00_01MemorySpace_00_01LayoutSpace_01_4.md)



[More...](#detailed-description)

* `#include <ipolar_poisson_like_solver.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef typename ddc::detail::convert\_type\_seq\_to\_discrete\_domain\_t&lt; [**batch\_tags**](classIPolarPoissonLikeSolver_3_01IdxRange_3_01ODims_8_8_8_01_4_00_01IdxRangeFull_00_01MemorySpace_00_01LayoutSpace_01_4.md#typedef-batch_tags) &gt; | [**batch\_idx\_range\_type**](#typedef-batch_idx_range_type)  <br>_The index range type describing the batch dimensions._  |
| typedef typename batch\_idx\_range\_type::discrete\_element\_type | [**batch\_index\_type**](#typedef-batch_index_type)  <br>_The index for indexing a batch dimension._  |
| typedef DConstField&lt; IdxRangeFull, MemorySpace, LayoutSpace &gt; | [**const\_field\_type**](#typedef-const_field_type)  <br>_The const Field type of the arguments to operator()._  |
| typedef DField&lt; IdxRangeFull, MemorySpace, LayoutSpace &gt; | [**field\_type**](#typedef-field_type)  <br>_The Field type of the arguments to operator()._  |
| typedef IdxRange&lt; ODims... &gt; | [**laplacian\_idx\_range\_type**](#typedef-laplacian_idx_range_type)  <br>_The type of the index range on which the equation is defined._  |
| typedef LayoutSpace | [**layout\_space**](#typedef-layout_space)  <br>_The layout space of the Fields passed to operator()._  |
| typedef MemorySpace | [**memory\_space**](#typedef-memory_space)  <br>_The space (CPU/GPU) where the Fields passed to operator() are saved._  |
| typedef std::conditional\_t&lt; ddc::type\_seq\_size\_v&lt; [**laplacian\_tags**](classIPolarPoissonLikeSolver_3_01IdxRange_3_01ODims_8_8_8_01_4_00_01IdxRangeFull_00_01MemorySpace_00_01LayoutSpace_01_4.md#typedef-laplacian_tags) &gt;==1, [**field\_type**](classIPolarPoissonLikeSolver_3_01IdxRange_3_01ODims_8_8_8_01_4_00_01IdxRangeFull_00_01MemorySpace_00_01LayoutSpace_01_4.md#typedef-field_type), [**VectorField**](classVectorField.md)&lt; double, IdxRangeFull, [**real\_laplacian\_tags**](classIPolarPoissonLikeSolver_3_01IdxRange_3_01ODims_8_8_8_01_4_00_01IdxRangeFull_00_01MemorySpace_00_01LayoutSpace_01_4.md#typedef-real_laplacian_tags), MemorySpace, LayoutSpace &gt; &gt; | [**vector\_field\_type**](#typedef-vector_field_type)  <br>_The type of the derivative of_ \(\phi\) _._ |




















## Public Functions

| Type | Name |
| ---: | :--- |
| virtual void | [**operator()**](#function-operator) ([**field\_type**](classIPolarPoissonLikeSolver_3_01IdxRange_3_01ODims_8_8_8_01_4_00_01IdxRangeFull_00_01MemorySpace_00_01LayoutSpace_01_4.md#typedef-field_type) phi, [**const\_field\_type**](classIPolarPoissonLikeSolver_3_01IdxRange_3_01ODims_8_8_8_01_4_00_01IdxRangeFull_00_01MemorySpace_00_01LayoutSpace_01_4.md#typedef-const_field_type) rho) const = 0<br>_An operator which calculates the solution_ \(\phi\) _to the Poisson-like equation: (1)_\(L\phi = - \nabla \cdot (\alpha \nabla \phi) + \beta \phi = \rho\) _, in_\(\Omega\) _,._ |
| virtual void | [**update\_coefficients**](#function-update_coefficients) ([**const\_field\_type**](classIPolarPoissonLikeSolver_3_01IdxRange_3_01ODims_8_8_8_01_4_00_01IdxRangeFull_00_01MemorySpace_00_01LayoutSpace_01_4.md#typedef-const_field_type) alpha, [**const\_field\_type**](classIPolarPoissonLikeSolver_3_01IdxRange_3_01ODims_8_8_8_01_4_00_01IdxRangeFull_00_01MemorySpace_00_01LayoutSpace_01_4.md#typedef-const_field_type) beta) = 0<br>_Update the coefficients_ \(alpha\) _and_\(beta\) _that define the equation._ |
| virtual  | [**~IPolarPoissonLikeSolver**](#function-ipolarpoissonlikesolver) () = default<br> |




## Protected Types

| Type | Name |
| ---: | :--- |
| typedef ddc::type\_seq\_remove\_t&lt; [**space\_tags**](classIPolarPoissonLikeSolver_3_01IdxRange_3_01ODims_8_8_8_01_4_00_01IdxRangeFull_00_01MemorySpace_00_01LayoutSpace_01_4.md#typedef-space_tags), [**laplacian\_tags**](classIPolarPoissonLikeSolver_3_01IdxRange_3_01ODims_8_8_8_01_4_00_01IdxRangeFull_00_01MemorySpace_00_01LayoutSpace_01_4.md#typedef-laplacian_tags) &gt; | [**batch\_tags**](#typedef-batch_tags)  <br>_The tags describing the batched dimensions._  |
| typedef ddc::detail::TypeSeq&lt; ODims... &gt; | [**laplacian\_tags**](#typedef-laplacian_tags)  <br>_The tags describing the discrete dimensions in the equation._  |
| typedef ddc::detail::TypeSeq&lt; typename ODims::continuous\_dimension\_type... &gt; | [**real\_laplacian\_tags**](#typedef-real_laplacian_tags)  <br>_The tags describing the real dimensions in the equation._  |
| typedef ddc::to\_type\_seq\_t&lt; IdxRangeFull &gt; | [**space\_tags**](#typedef-space_tags)  <br>_The tags describing the dimensions of the index range on which the operator acts._  |






## Protected Static Attributes

| Type | Name |
| ---: | :--- |
|  constexpr bool | [**using\_vector\_field**](#variable-using_vector_field)   = `ddc::type\_seq\_size\_v&lt;[**laplacian\_tags**](classIPolarPoissonLikeSolver_3_01IdxRange_3_01ODims_8_8_8_01_4_00_01IdxRangeFull_00_01MemorySpace_00_01LayoutSpace_01_4.md#typedef-laplacian_tags)&gt; == 1`<br>_Indicates whether the gradient is represented by a_ [_**VectorField**_](classVectorField.md) _or a Field._ |


















## Detailed Description


An abstract class from which a Poisson-like solver can inherit. Classes inheriting from this must implement a way to solve the following equation: (1) \(L\phi = - \nabla \cdot (\alpha \nabla \phi) + \beta \phi = \rho\), in \(\Omega\),


\(\phi = 0\), on \(\partial \Omega\),




**Template parameters:**


* `IdxRangeLaplacian` The index range on which the equation is defined. 
* `IdxRangeFull` The index range on which the operator() acts. This is equal to the IdxRangeLaplacian plus any batched dimensions. 
* `MemorySpace` The space (CPU/GPU) where the Fields passed to operator() are saved. 
* `LayoutSpace` The layout space of the Fields passed to operator(). 




    
## Public Types Documentation




### typedef batch\_idx\_range\_type 

_The index range type describing the batch dimensions._ 
```C++
using IPolarPoissonLikeSolver< IdxRange< ODims... >, IdxRangeFull, MemorySpace, LayoutSpace >::batch_idx_range_type =  typename ddc::detail::convert_type_seq_to_discrete_domain_t<batch_tags>;
```




<hr>



### typedef batch\_index\_type 

_The index for indexing a batch dimension._ 
```C++
using IPolarPoissonLikeSolver< IdxRange< ODims... >, IdxRangeFull, MemorySpace, LayoutSpace >::batch_index_type =  typename batch_idx_range_type::discrete_element_type;
```




<hr>



### typedef const\_field\_type 

_The const Field type of the arguments to operator()._ 
```C++
using IPolarPoissonLikeSolver< IdxRange< ODims... >, IdxRangeFull, MemorySpace, LayoutSpace >::const_field_type =  DConstField<IdxRangeFull, MemorySpace, LayoutSpace>;
```




<hr>



### typedef field\_type 

_The Field type of the arguments to operator()._ 
```C++
using IPolarPoissonLikeSolver< IdxRange< ODims... >, IdxRangeFull, MemorySpace, LayoutSpace >::field_type =  DField<IdxRangeFull, MemorySpace, LayoutSpace>;
```




<hr>



### typedef laplacian\_idx\_range\_type 

_The type of the index range on which the equation is defined._ 
```C++
using IPolarPoissonLikeSolver< IdxRange< ODims... >, IdxRangeFull, MemorySpace, LayoutSpace >::laplacian_idx_range_type =  IdxRange<ODims...>;
```




<hr>



### typedef layout\_space 

_The layout space of the Fields passed to operator()._ 
```C++
using IPolarPoissonLikeSolver< IdxRange< ODims... >, IdxRangeFull, MemorySpace, LayoutSpace >::layout_space =  LayoutSpace;
```




<hr>



### typedef memory\_space 

_The space (CPU/GPU) where the Fields passed to operator() are saved._ 
```C++
using IPolarPoissonLikeSolver< IdxRange< ODims... >, IdxRangeFull, MemorySpace, LayoutSpace >::memory_space =  MemorySpace;
```




<hr>



### typedef vector\_field\_type 

_The type of the derivative of_ \(\phi\) _._
```C++
using IPolarPoissonLikeSolver< IdxRange< ODims... >, IdxRangeFull, MemorySpace, LayoutSpace >::vector_field_type =  std::conditional_t< ddc::type_seq_size_v<laplacian_tags> == 1, field_type, VectorField<double, IdxRangeFull, real_laplacian_tags, MemorySpace, LayoutSpace> >;
```




<hr>
## Public Functions Documentation




### function operator() 

_An operator which calculates the solution_ \(\phi\) _to the Poisson-like equation: (1)_\(L\phi = - \nabla \cdot (\alpha \nabla \phi) + \beta \phi = \rho\) _, in_\(\Omega\) _,._
```C++
virtual void IPolarPoissonLikeSolver< IdxRange< ODims... >, IdxRangeFull, MemorySpace, LayoutSpace >::operator() (
    field_type phi,
    const_field_type rho
) const = 0
```



\(\phi = 0\), on \(\partial \Omega\),




**Parameters:**


* `phi` The solution to Poisson's equation. 
* `rho` The right-hand side of Poisson's equation. 




        

<hr>



### function update\_coefficients 

_Update the coefficients_ \(alpha\) _and_\(beta\) _that define the equation._
```C++
virtual void IPolarPoissonLikeSolver< IdxRange< ODims... >, IdxRangeFull, MemorySpace, LayoutSpace >::update_coefficients (
    const_field_type alpha,
    const_field_type beta
) = 0
```





**Parameters:**


* `alpha` The values of alpha at the grid points. 
* `beta` The values of beta at the grid points. 




        

<hr>



### function ~IPolarPoissonLikeSolver 

```C++
virtual IPolarPoissonLikeSolver< IdxRange< ODims... >, IdxRangeFull, MemorySpace, LayoutSpace >::~IPolarPoissonLikeSolver () = default
```




<hr>
## Protected Types Documentation




### typedef batch\_tags 

_The tags describing the batched dimensions._ 
```C++
using IPolarPoissonLikeSolver< IdxRange< ODims... >, IdxRangeFull, MemorySpace, LayoutSpace >::batch_tags =  ddc::type_seq_remove_t<space_tags, laplacian_tags>;
```




<hr>



### typedef laplacian\_tags 

_The tags describing the discrete dimensions in the equation._ 
```C++
using IPolarPoissonLikeSolver< IdxRange< ODims... >, IdxRangeFull, MemorySpace, LayoutSpace >::laplacian_tags =  ddc::detail::TypeSeq<ODims...>;
```




<hr>



### typedef real\_laplacian\_tags 

_The tags describing the real dimensions in the equation._ 
```C++
using IPolarPoissonLikeSolver< IdxRange< ODims... >, IdxRangeFull, MemorySpace, LayoutSpace >::real_laplacian_tags =  ddc::detail::TypeSeq<typename ODims::continuous_dimension_type...>;
```




<hr>



### typedef space\_tags 

_The tags describing the dimensions of the index range on which the operator acts._ 
```C++
using IPolarPoissonLikeSolver< IdxRange< ODims... >, IdxRangeFull, MemorySpace, LayoutSpace >::space_tags =  ddc::to_type_seq_t<IdxRangeFull>;
```




<hr>
## Protected Static Attributes Documentation




### variable using\_vector\_field 

_Indicates whether the gradient is represented by a_ [_**VectorField**_](classVectorField.md) _or a Field._
```C++
constexpr bool IPolarPoissonLikeSolver< IdxRange< ODims... >, IdxRangeFull, MemorySpace, LayoutSpace >::using_vector_field;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/pde_solvers/ipolar_poisson_like_solver.hpp`

