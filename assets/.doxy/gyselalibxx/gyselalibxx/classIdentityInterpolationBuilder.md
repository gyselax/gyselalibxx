

# Class IdentityInterpolationBuilder

**template &lt;class ExecSpace, class MemorySpace, class InterpolationDDim, class Basis&gt;**



[**ClassList**](annotated.md) **>** [**IdentityInterpolationBuilder**](classIdentityInterpolationBuilder.md)



_A builder class for copying data._ [More...](#detailed-description)

* `#include <identity_interpolation_builder.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef typename Basis::template Impl&lt; Basis, MemorySpace &gt;::knot\_grid | [**basis\_domain\_type**](#typedef-basis_domain_type)  <br>_The grid on which the interpolation coefficients should be provided._  |
| typedef ddc::remove\_dims\_of\_t&lt; BatchedInterpolationGrid, [**interpolation\_discrete\_dimension\_type**](classIdentityInterpolationBuilder.md#typedef-interpolation_discrete_dimension_type) &gt; | [**batch\_domain\_type**](#typedef-batch_domain_type)  <br>_The type of the batch domain (obtained by removing the dimension of interest from the whole domain)._  |
| typedef ddc::replace\_dim\_of\_t&lt; BatchedInterpolationGrid, [**interpolation\_discrete\_dimension\_type**](classIdentityInterpolationBuilder.md#typedef-interpolation_discrete_dimension_type), [**basis\_domain\_type**](classIdentityInterpolationBuilder.md#typedef-basis_domain_type) &gt; | [**batched\_basis\_domain\_type**](#typedef-batched_basis_domain_type)  <br>_The type of the whole interpolation domain (cartesian product of 1D interpolation domain and batch domain) preserving the underlying memory layout (order of dimensions)._  |
| typedef ddc::remove\_dims\_of\_t&lt; BatchedInterpolationGrid, [**interpolation\_discrete\_dimension\_type**](classIdentityInterpolationBuilder.md#typedef-interpolation_discrete_dimension_type) &gt; | [**batched\_derivs\_domain\_type**](#typedef-batched_derivs_domain_type)  <br>_The typeof the derivatives._  |
| typedef BatchedInterpolationGrid | [**batched\_interpolation\_domain\_type**](#typedef-batched_interpolation_domain_type)  <br>_The type of the whole domain representing interpolation points._  |
| typedef typename InterpolationDDim::continuous\_dimension\_type | [**continuous\_dimension\_type**](#typedef-continuous_dimension_type)  <br>_The type of the interpolation continuous dimension (continuous dimension of interest) used by this class._  |
| typedef ExecSpace | [**exec\_space**](#typedef-exec_space)  <br>_The type of the Kokkos execution space used by this class._  |
| typedef InterpolationDDim | [**interpolation\_discrete\_dimension\_type**](#typedef-interpolation_discrete_dimension_type)  <br>_The type of the interpolation discrete dimension (discrete dimension of interest) used by this class._  |
| typedef IdxRange&lt; [**interpolation\_discrete\_dimension\_type**](classIdentityInterpolationBuilder.md#typedef-interpolation_discrete_dimension_type) &gt; | [**interpolation\_domain\_type**](#typedef-interpolation_domain_type)  <br>_The type of the domain for the 1D interpolation mesh used by this class._  |
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
|  void | [**operator()**](#function-operator) (Field&lt; DataType, [**batched\_basis\_domain\_type**](classIdentityInterpolationBuilder.md#typedef-batched_basis_domain_type)&lt; BatchedInterpolationGrid &gt;, [**memory\_space**](classIdentityInterpolationBuilder.md#typedef-memory_space), Layout &gt; coeffs, ConstField&lt; DataType, BatchedInterpolationGrid, [**memory\_space**](classIdentityInterpolationBuilder.md#typedef-memory_space), Layout &gt; vals, std::optional&lt; ConstField&lt; DataType, [**batched\_derivs\_domain\_type**](classIdentityInterpolationBuilder.md#typedef-batched_derivs_domain_type)&lt; BatchedInterpolationGrid &gt;, [**memory\_space**](classIdentityInterpolationBuilder.md#typedef-memory_space), Layout &gt; &gt; derivs\_xmin=std::nullopt, std::optional&lt; ConstField&lt; DataType, [**batched\_derivs\_domain\_type**](classIdentityInterpolationBuilder.md#typedef-batched_derivs_domain_type)&lt; BatchedInterpolationGrid &gt;, [**memory\_space**](classIdentityInterpolationBuilder.md#typedef-memory_space), Layout &gt; &gt; derivs\_xmax=std::nullopt) const<br>_Compute the interpolation coefficients for a function._  |




























## Detailed Description


A class which contains an operator () which can be used to build an interpolation of a function. This class handles the case where no calculations are necessary and the data simply needs to be copied.




**Template parameters:**


* `ExecSpace` The Kokkos execution space on which the spline approximation is performed. 
* `MemorySpace` The Kokkos memory space on which the data (interpolation function and splines coefficients) is stored. 
* `InterpolationDDim` The discrete dimension on which interpolation points are defined. 
* `Basis` The basis on which the interpolation is constructed. 




    
## Public Types Documentation




### typedef basis\_domain\_type 

_The grid on which the interpolation coefficients should be provided._ 
```C++
using IdentityInterpolationBuilder< ExecSpace, MemorySpace, InterpolationDDim, Basis >::basis_domain_type =  typename Basis::template Impl<Basis, MemorySpace>::knot_grid;
```




<hr>



### typedef batch\_domain\_type 

_The type of the batch domain (obtained by removing the dimension of interest from the whole domain)._ 
```C++
using IdentityInterpolationBuilder< ExecSpace, MemorySpace, InterpolationDDim, Basis >::batch_domain_type =  ddc:: remove_dims_of_t<BatchedInterpolationGrid, interpolation_discrete_dimension_type>;
```





**Template parameters:**


* `The` batched discrete domain on which the interpolation points are defined.

Example: For batched\_interpolation\_domain\_type = DiscreteDomain&lt;X,Y,Z&gt; and a dimension of interest [**Y**](structY.md), this is DiscreteDomain&lt;X,Z&gt; 


        

<hr>



### typedef batched\_basis\_domain\_type 

_The type of the whole interpolation domain (cartesian product of 1D interpolation domain and batch domain) preserving the underlying memory layout (order of dimensions)._ 
```C++
using IdentityInterpolationBuilder< ExecSpace, MemorySpace, InterpolationDDim, Basis >::batched_basis_domain_type =  ddc::replace_dim_of_t< BatchedInterpolationGrid, interpolation_discrete_dimension_type, basis_domain_type>;
```





**Template parameters:**


* `The` batched discrete domain on which the interpolation points are defined.

Example: For batched\_interpolation\_domain\_type = DiscreteDomain&lt;X,Y,Z&gt; and a dimension of interest [**Y**](structY.md) (associated to a Basis tag LagrangeBasisY), this is DiscreteDomain&lt;X,LagrangeBasisY,Z&gt;. 


        

<hr>



### typedef batched\_derivs\_domain\_type 

_The typeof the derivatives._ 
```C++
using IdentityInterpolationBuilder< ExecSpace, MemorySpace, InterpolationDDim, Basis >::batched_derivs_domain_type =  ddc:: remove_dims_of_t<BatchedInterpolationGrid, interpolation_discrete_dimension_type>;
```



The type of the derivatives that need to be provided to this method. No derivatives are required but this is included for interoperability with other interpolation builder classes.




**Template parameters:**


* `The` batched discrete domain on which the interpolation points are defined. 




        

<hr>



### typedef batched\_interpolation\_domain\_type 

_The type of the whole domain representing interpolation points._ 
```C++
using IdentityInterpolationBuilder< ExecSpace, MemorySpace, InterpolationDDim, Basis >::batched_interpolation_domain_type =  BatchedInterpolationGrid;
```





**Template parameters:**


* `The` batched discrete domain on which the interpolation points are defined. 




        

<hr>



### typedef continuous\_dimension\_type 

_The type of the interpolation continuous dimension (continuous dimension of interest) used by this class._ 
```C++
using IdentityInterpolationBuilder< ExecSpace, MemorySpace, InterpolationDDim, Basis >::continuous_dimension_type =  typename InterpolationDDim::continuous_dimension_type;
```




<hr>



### typedef exec\_space 

_The type of the Kokkos execution space used by this class._ 
```C++
using IdentityInterpolationBuilder< ExecSpace, MemorySpace, InterpolationDDim, Basis >::exec_space =  ExecSpace;
```




<hr>



### typedef interpolation\_discrete\_dimension\_type 

_The type of the interpolation discrete dimension (discrete dimension of interest) used by this class._ 
```C++
using IdentityInterpolationBuilder< ExecSpace, MemorySpace, InterpolationDDim, Basis >::interpolation_discrete_dimension_type =  InterpolationDDim;
```




<hr>



### typedef interpolation\_domain\_type 

_The type of the domain for the 1D interpolation mesh used by this class._ 
```C++
using IdentityInterpolationBuilder< ExecSpace, MemorySpace, InterpolationDDim, Basis >::interpolation_domain_type =  IdxRange<interpolation_discrete_dimension_type>;
```




<hr>



### typedef memory\_space 

_The type of the Kokkos memory space used by this class._ 
```C++
using IdentityInterpolationBuilder< ExecSpace, MemorySpace, InterpolationDDim, Basis >::memory_space =  MemorySpace;
```




<hr>
## Public Static Attributes Documentation




### variable s\_nbc\_xmax 

_The number of equations defining the boundary condition at the upper bound._ 
```C++
constexpr int IdentityInterpolationBuilder< ExecSpace, MemorySpace, InterpolationDDim, Basis >::s_nbc_xmax;
```




<hr>



### variable s\_nbc\_xmin 

_The number of equations defining the boundary condition at the lower bound._ 
```C++
constexpr int IdentityInterpolationBuilder< ExecSpace, MemorySpace, InterpolationDDim, Basis >::s_nbc_xmin;
```




<hr>
## Public Functions Documentation




### function IdentityInterpolationBuilder 

```C++
IdentityInterpolationBuilder::IdentityInterpolationBuilder () = default
```




<hr>



### function operator() 

_Compute the interpolation coefficients for a function._ 
```C++
template<class DataType, class Layout, class BatchedInterpolationGrid>
inline void IdentityInterpolationBuilder::operator() (
    Field< DataType, batched_basis_domain_type < BatchedInterpolationGrid >, memory_space , Layout > coeffs,
    ConstField< DataType, BatchedInterpolationGrid, memory_space , Layout > vals,
    std::optional< ConstField< DataType, batched_derivs_domain_type < BatchedInterpolationGrid >, memory_space , Layout > > derivs_xmin=std::nullopt,
    std::optional< ConstField< DataType, batched_derivs_domain_type < BatchedInterpolationGrid >, memory_space , Layout > > derivs_xmax=std::nullopt
) const
```



No calculations are necessary and the data simply needs to be copied from vals to coeffs.




**Parameters:**


* `coeffs` The coefficients of the spline computed by this SplineBuilder. 
* `vals` The values of the function on the interpolation mesh. 
* `derivs_xmin` The values of the derivatives at the lower boundary (unused in this class). 
* `derivs_xmax` The values of the derivatives at the upper boundary (unused in this class). 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/interpolation/identity_interpolation_builder.hpp`

