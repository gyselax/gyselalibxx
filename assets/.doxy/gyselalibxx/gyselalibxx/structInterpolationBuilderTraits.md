

# Struct InterpolationBuilderTraits

**template &lt;class Builder&gt;**



[**ClassList**](annotated.md) **>** [**InterpolationBuilderTraits**](structInterpolationBuilderTraits.md)



_A traits struct for accessing type aliases of an interpolation builder._ [More...](#detailed-description)

* `#include <i_interpolation_builder.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef typename Builder::template [**batched\_basis\_idx\_range\_type**](structInterpolationBuilderTraits.md#typedef-batched_basis_idx_range_type)&lt; IdxRangeBatchedInterpolation &gt; | [**batched\_basis\_idx\_range\_type**](#typedef-batched_basis_idx_range_type)  <br>_Batched domain with the interpolation grid(s) replaced by the basis grid(s)._  |
| typedef typename Builder::template [**batched\_derivs\_idx\_range\_type**](structInterpolationBuilderTraits.md#typedef-batched_derivs_idx_range_type)&lt; IdxRangeBatchedInterpolation &gt; | [**batched\_derivs\_idx\_range\_type**](#typedef-batched_derivs_idx_range_type)  <br>_Batched domain with the interpolation grid replaced by the deriv type._  |
| typedef typename Builder::coeff\_idx\_range\_type | [**coeff\_idx\_range\_type**](#typedef-coeff_idx_range_type)  <br>_The index range for the interpolation coefficients._  |
| typedef typename Builder::data\_type | [**data\_type**](#typedef-data_type)  <br>_The data type that the data is saved on._  |
| typedef typename Builder::interpolation\_idx\_range\_type | [**interpolation\_idx\_range\_type**](#typedef-interpolation_idx_range_type)  <br>_The ND index range for the interpolation mesh._  |






















## Public Static Functions

| Type | Name |
| ---: | :--- |
|  constexpr std::size\_t | [**rank**](#function-rank) () <br>_The number of interpolation dimensions._  |


























## Detailed Description


The primary template delegates to the builder's own type aliases, so any class that defines them directly (e.g. [**IdentityInterpolationBuilder**](classIdentityInterpolationBuilder.md)) satisfies the InterpolationBuilder concept without specialisation.


Specialise this struct to adapt external builders whose alias names differ from the convention (e.g. ddc::SplineBuilder).


Defines: Type aliases:
* data\_type
* interpolation\_idx\_range\_type
* coeff\_idx\_range\_type Static functions:
* [**rank()**](structInterpolationBuilderTraits.md#function-rank) Type calculators:
* batched\_basis\_idx\_range\_type
* batched\_derivs\_idx\_range\_type (1D builders only; not required by the concept)






**Template parameters:**


* `Builder` The interpolation builder type. 




    
## Public Types Documentation




### typedef batched\_basis\_idx\_range\_type 

_Batched domain with the interpolation grid(s) replaced by the basis grid(s)._ 
```C++
using InterpolationBuilderTraits< Builder >::batched_basis_idx_range_type =  typename Builder::template batched_basis_idx_range_type<IdxRangeBatchedInterpolation>;
```




<hr>



### typedef batched\_derivs\_idx\_range\_type 

_Batched domain with the interpolation grid replaced by the deriv type._ 
```C++
using InterpolationBuilderTraits< Builder >::batched_derivs_idx_range_type =  typename Builder::template batched_derivs_idx_range_type<IdxRangeBatchedInterpolation>;
```




<hr>



### typedef coeff\_idx\_range\_type 

_The index range for the interpolation coefficients._ 
```C++
using InterpolationBuilderTraits< Builder >::coeff_idx_range_type =  typename Builder::coeff_idx_range_type;
```




<hr>



### typedef data\_type 

_The data type that the data is saved on._ 
```C++
using InterpolationBuilderTraits< Builder >::data_type =  typename Builder::data_type;
```




<hr>



### typedef interpolation\_idx\_range\_type 

_The ND index range for the interpolation mesh._ 
```C++
using InterpolationBuilderTraits< Builder >::interpolation_idx_range_type =  typename Builder::interpolation_idx_range_type;
```




<hr>
## Public Static Functions Documentation




### function rank 

_The number of interpolation dimensions._ 
```C++
static inline constexpr std::size_t InterpolationBuilderTraits::rank () 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/interpolation/i_interpolation_builder.hpp`

