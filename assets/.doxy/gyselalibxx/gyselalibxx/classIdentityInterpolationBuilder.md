

# Class IdentityInterpolationBuilder

**template &lt;class ExecSpace, class MemorySpace, class DataType, class InterpolationGrid, class Basis&gt;**



[**ClassList**](annotated.md) **>** [**IdentityInterpolationBuilder**](classIdentityInterpolationBuilder.md)



_A builder class for copying data._ [More...](#detailed-description)

* `#include <identity_interpolation_builder.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef typename Basis::template Impl&lt; Basis, MemorySpace &gt;::knot\_grid | [**basis\_domain\_type**](#typedef-basis_domain_type)  <br>_The grid on which the interpolation coefficients should be provided._  |
| typedef ddc::replace\_dim\_of\_t&lt; BatchedInterpolationIdxRange, [**interpolation\_grid\_type**](classIdentityInterpolationBuilder.md#typedef-interpolation_grid_type), [**basis\_domain\_type**](classIdentityInterpolationBuilder.md#typedef-basis_domain_type) &gt; | [**batched\_basis\_idx\_range\_type**](#typedef-batched_basis_idx_range_type)  <br>_The batched domain type with interpolation\_grid\_type replaced by basis\_domain\_type._  |
| typedef ddc::replace\_dim\_of\_t&lt; BatchedInterpolationIdxRange, [**interpolation\_grid\_type**](classIdentityInterpolationBuilder.md#typedef-interpolation_grid_type), [**deriv\_type**](classIdentityInterpolationBuilder.md#typedef-deriv_type) &gt; | [**batched\_derivs\_idx\_range\_type**](#typedef-batched_derivs_idx_range_type)  <br>_The batched domain type with interpolation\_grid\_type replaced by deriv\_type._  |
| typedef typename InterpolationGrid::continuous\_dimension\_type | [**continuous\_dimension\_type**](#typedef-continuous_dimension_type)  <br>_The type of the interpolation continuous dimension (continuous dimension of interest) used by this class._  |
| typedef DataType | [**data\_type**](#typedef-data_type)  <br>_The data type that the data is saved on._  |
| typedef ddc::Deriv&lt; [**continuous\_dimension\_type**](classIdentityInterpolationBuilder.md#typedef-continuous_dimension_type) &gt; | [**deriv\_type**](#typedef-deriv_type)  <br>_The type of the Deriv dimension at the boundaries._  |
| typedef ExecSpace | [**exec\_space**](#typedef-exec_space)  <br>_The type of the Kokkos execution space used by this class._  |
| typedef InterpolationGrid | [**interpolation\_grid\_type**](#typedef-interpolation_grid_type)  <br>_The type of the interpolation discrete dimension (discrete dimension of interest) used by this class._  |
| typedef IdxRange&lt; [**interpolation\_grid\_type**](classIdentityInterpolationBuilder.md#typedef-interpolation_grid_type) &gt; | [**interpolation\_idx\_range\_type**](#typedef-interpolation_idx_range_type)  <br>_The type of the domain for the 1D interpolation mesh used by this class._  |
| typedef MemorySpace | [**memory\_space**](#typedef-memory_space)  <br>_The type of the Kokkos memory space used by this class._  |






## Public Static Attributes

| Type | Name |
| ---: | :--- |
|  constexpr int | [**s\_nbc\_xmax**](#variable-s_nbc_xmax)   = `0`<br>_The number of equations defining the boundary condition at the upper bound._  |
|  constexpr int | [**s\_nbc\_xmin**](#variable-s_nbc_xmin)   = `0`<br>_The number of equations defining the boundary condition at the lower bound._  |














## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**IdentityInterpolationBuilder**](#function-identityinterpolationbuilder) () = default<br> |
|  [**batched\_basis\_idx\_range\_type**](classIdentityInterpolationBuilder.md#typedef-batched_basis_idx_range_type)&lt; BatchedInterpolationIdxRange &gt; | [**batched\_basis\_idx\_range**](#function-batched_basis_idx_range) (BatchedInterpolationIdxRange const & batched\_interpolation\_domain) noexcept const<br>_Get the whole domain on which the interpolation coefficients are defined._  |
|  [**batched\_derivs\_idx\_range\_type**](classIdentityInterpolationBuilder.md#typedef-batched_derivs_idx_range_type)&lt; BatchedInterpolationIdxRange &gt; | [**batched\_derivs\_xmax\_domain**](#function-batched_derivs_xmax_domain) (BatchedInterpolationIdxRange const & batched\_interpolation\_domain) noexcept const<br>_Get the whole domain on which derivatives on upper boundary are defined._  |
|  [**batched\_derivs\_idx\_range\_type**](classIdentityInterpolationBuilder.md#typedef-batched_derivs_idx_range_type)&lt; BatchedInterpolationIdxRange &gt; | [**batched\_derivs\_xmin\_domain**](#function-batched_derivs_xmin_domain) (BatchedInterpolationIdxRange const & batched\_interpolation\_domain) noexcept const<br>_Get the whole domain on which derivatives on lower boundary are defined._  |
|  void | [**operator()**](#function-operator) (Field&lt; DataType, [**batched\_basis\_idx\_range\_type**](classIdentityInterpolationBuilder.md#typedef-batched_basis_idx_range_type)&lt; BatchedInterpolationIdxRange &gt;, [**memory\_space**](classIdentityInterpolationBuilder.md#typedef-memory_space) &gt; coeffs, ConstField&lt; DataType, BatchedInterpolationIdxRange, [**memory\_space**](classIdentityInterpolationBuilder.md#typedef-memory_space) &gt; vals, std::optional&lt; ConstField&lt; DataType, [**batched\_derivs\_idx\_range\_type**](classIdentityInterpolationBuilder.md#typedef-batched_derivs_idx_range_type)&lt; BatchedInterpolationIdxRange &gt;, [**memory\_space**](classIdentityInterpolationBuilder.md#typedef-memory_space) &gt; &gt; derivs\_xmin=std::nullopt, std::optional&lt; ConstField&lt; DataType, [**batched\_derivs\_idx\_range\_type**](classIdentityInterpolationBuilder.md#typedef-batched_derivs_idx_range_type)&lt; BatchedInterpolationIdxRange &gt;, [**memory\_space**](classIdentityInterpolationBuilder.md#typedef-memory_space) &gt; &gt; derivs\_xmax=std::nullopt) const<br>_Compute the interpolation coefficients for a function._  |




























## Detailed Description


A class which contains an operator () which can be used to build an interpolation of a function. This class handles the case where no calculations are necessary and the data simply needs to be copied.




**Template parameters:**


* `ExecSpace` The Kokkos execution space on which the spline approximation is performed. 
* `MemorySpace` The Kokkos memory space on which the data (interpolation function and splines coefficients) is stored. 
* `DataType` The data type of the field values and coefficients. 
* `InterpolationGrid` The discrete dimension on which interpolation points are defined. 
* `Basis` The basis on which the interpolation is constructed. 




    
## Public Types Documentation




### typedef basis\_domain\_type 

_The grid on which the interpolation coefficients should be provided._ 
```C++
using IdentityInterpolationBuilder< ExecSpace, MemorySpace, DataType, InterpolationGrid, Basis >::basis_domain_type =  typename Basis::template Impl<Basis, MemorySpace>::knot_grid;
```




<hr>



### typedef batched\_basis\_idx\_range\_type 

_The batched domain type with interpolation\_grid\_type replaced by basis\_domain\_type._ 
```C++
using IdentityInterpolationBuilder< ExecSpace, MemorySpace, DataType, InterpolationGrid, Basis >::batched_basis_idx_range_type =  ddc::replace_dim_of_t< BatchedInterpolationIdxRange, interpolation_grid_type, basis_domain_type>;
```




<hr>



### typedef batched\_derivs\_idx\_range\_type 

_The batched domain type with interpolation\_grid\_type replaced by deriv\_type._ 
```C++
using IdentityInterpolationBuilder< ExecSpace, MemorySpace, DataType, InterpolationGrid, Basis >::batched_derivs_idx_range_type =  ddc:: replace_dim_of_t<BatchedInterpolationIdxRange, interpolation_grid_type, deriv_type>;
```




<hr>



### typedef continuous\_dimension\_type 

_The type of the interpolation continuous dimension (continuous dimension of interest) used by this class._ 
```C++
using IdentityInterpolationBuilder< ExecSpace, MemorySpace, DataType, InterpolationGrid, Basis >::continuous_dimension_type =  typename InterpolationGrid::continuous_dimension_type;
```




<hr>



### typedef data\_type 

_The data type that the data is saved on._ 
```C++
using IdentityInterpolationBuilder< ExecSpace, MemorySpace, DataType, InterpolationGrid, Basis >::data_type =  DataType;
```




<hr>



### typedef deriv\_type 

_The type of the Deriv dimension at the boundaries._ 
```C++
using IdentityInterpolationBuilder< ExecSpace, MemorySpace, DataType, InterpolationGrid, Basis >::deriv_type =  ddc::Deriv<continuous_dimension_type>;
```




<hr>



### typedef exec\_space 

_The type of the Kokkos execution space used by this class._ 
```C++
using IdentityInterpolationBuilder< ExecSpace, MemorySpace, DataType, InterpolationGrid, Basis >::exec_space =  ExecSpace;
```




<hr>



### typedef interpolation\_grid\_type 

_The type of the interpolation discrete dimension (discrete dimension of interest) used by this class._ 
```C++
using IdentityInterpolationBuilder< ExecSpace, MemorySpace, DataType, InterpolationGrid, Basis >::interpolation_grid_type =  InterpolationGrid;
```




<hr>



### typedef interpolation\_idx\_range\_type 

_The type of the domain for the 1D interpolation mesh used by this class._ 
```C++
using IdentityInterpolationBuilder< ExecSpace, MemorySpace, DataType, InterpolationGrid, Basis >::interpolation_idx_range_type =  IdxRange<interpolation_grid_type>;
```




<hr>



### typedef memory\_space 

_The type of the Kokkos memory space used by this class._ 
```C++
using IdentityInterpolationBuilder< ExecSpace, MemorySpace, DataType, InterpolationGrid, Basis >::memory_space =  MemorySpace;
```




<hr>
## Public Static Attributes Documentation




### variable s\_nbc\_xmax 

_The number of equations defining the boundary condition at the upper bound._ 
```C++
constexpr int IdentityInterpolationBuilder< ExecSpace, MemorySpace, DataType, InterpolationGrid, Basis >::s_nbc_xmax;
```




<hr>



### variable s\_nbc\_xmin 

_The number of equations defining the boundary condition at the lower bound._ 
```C++
constexpr int IdentityInterpolationBuilder< ExecSpace, MemorySpace, DataType, InterpolationGrid, Basis >::s_nbc_xmin;
```




<hr>
## Public Functions Documentation




### function IdentityInterpolationBuilder 

```C++
IdentityInterpolationBuilder::IdentityInterpolationBuilder () = default
```




<hr>



### function batched\_basis\_idx\_range 

_Get the whole domain on which the interpolation coefficients are defined._ 
```C++
template<class BatchedInterpolationIdxRange>
inline batched_basis_idx_range_type < BatchedInterpolationIdxRange > IdentityInterpolationBuilder::batched_basis_idx_range (
    BatchedInterpolationIdxRange const & batched_interpolation_domain
) noexcept const
```





**Parameters:**


* `batched_interpolation_domain` The whole domain on which the interpolation points are defined.



**Returns:**

The domain for the interpolation coefficients. 





        

<hr>



### function batched\_derivs\_xmax\_domain 

_Get the whole domain on which derivatives on upper boundary are defined._ 
```C++
template<class BatchedInterpolationIdxRange>
inline batched_derivs_idx_range_type < BatchedInterpolationIdxRange > IdentityInterpolationBuilder::batched_derivs_xmax_domain (
    BatchedInterpolationIdxRange const & batched_interpolation_domain
) noexcept const
```



This is only used with BoundCond::HERMITE boundary conditions.




**Parameters:**


* `batched_interpolation_domain` The whole domain on which the interpolation points are defined.



**Returns:**

The domain for the Derivs values. 





        

<hr>



### function batched\_derivs\_xmin\_domain 

_Get the whole domain on which derivatives on lower boundary are defined._ 
```C++
template<class BatchedInterpolationIdxRange>
inline batched_derivs_idx_range_type < BatchedInterpolationIdxRange > IdentityInterpolationBuilder::batched_derivs_xmin_domain (
    BatchedInterpolationIdxRange const & batched_interpolation_domain
) noexcept const
```



This is only used with BoundCond::HERMITE boundary conditions.




**Parameters:**


* `batched_interpolation_domain` The whole domain on which the interpolation points are defined.



**Returns:**

The domain for the Derivs values. 





        

<hr>



### function operator() 

_Compute the interpolation coefficients for a function._ 
```C++
template<class BatchedInterpolationIdxRange>
inline void IdentityInterpolationBuilder::operator() (
    Field< DataType, batched_basis_idx_range_type < BatchedInterpolationIdxRange >, memory_space > coeffs,
    ConstField< DataType, BatchedInterpolationIdxRange, memory_space > vals,
    std::optional< ConstField< DataType, batched_derivs_idx_range_type < BatchedInterpolationIdxRange >, memory_space > > derivs_xmin=std::nullopt,
    std::optional< ConstField< DataType, batched_derivs_idx_range_type < BatchedInterpolationIdxRange >, memory_space > > derivs_xmax=std::nullopt
) const
```



No calculations are necessary and the data simply needs to be copied from vals to coeffs.




**Parameters:**


* `coeffs` The coefficients of the interpolation computed by this builder. 
* `vals` The values of the function on the interpolation mesh. 
* `derivs_xmin` The values of the derivatives at the lower boundary (unused in this class). 
* `derivs_xmax` The values of the derivatives at the upper boundary (unused in this class). 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/interpolation/identity_interpolation_builder.hpp`

