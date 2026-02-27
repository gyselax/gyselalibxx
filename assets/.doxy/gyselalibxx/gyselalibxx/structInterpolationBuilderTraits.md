

# Struct InterpolationBuilderTraits

**template &lt;class Builder&gt;**



[**ClassList**](annotated.md) **>** [**InterpolationBuilderTraits**](structInterpolationBuilderTraits.md)



_A traits struct for accessing type aliases of an interpolation builder._ [More...](#detailed-description)

* `#include <i_interpolation_builder.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef typename Builder::basis\_domain\_type | [**basis\_domain\_type**](#typedef-basis_domain_type)  <br>_The discrete dimension for the interpolation coefficients._  |
| typedef typename Builder::template [**batched\_basis\_idx\_range\_type**](structInterpolationBuilderTraits.md#typedef-batched_basis_idx_range_type)&lt; IdxRangeBatchedInterpolation &gt; | [**batched\_basis\_idx\_range\_type**](#typedef-batched_basis_idx_range_type)  <br>_Batched domain with interpolation\_grid\_type replaced by basis\_domain\_type._  |
| typedef typename Builder::template [**batched\_derivs\_idx\_range\_type**](structInterpolationBuilderTraits.md#typedef-batched_derivs_idx_range_type)&lt; IdxRangeBatchedInterpolation &gt; | [**batched\_derivs\_idx\_range\_type**](#typedef-batched_derivs_idx_range_type)  <br>_Batched domain with interpolation\_grid\_type replaced by deriv\_type._  |
| typedef typename Builder::data\_type | [**data\_type**](#typedef-data_type)  <br>_The data type that the data is saved on._  |
| typedef typename Builder::interpolation\_grid\_type | [**interpolation\_grid\_type**](#typedef-interpolation_grid_type)  <br>_The discrete grid on which interpolation values are given._  |
| typedef typename Builder::interpolation\_idx\_range\_type | [**interpolation\_idx\_range\_type**](#typedef-interpolation_idx_range_type)  <br>_The 1D index range for the interpolation mesh._  |
















































## Detailed Description


The primary template delegates to the builder's own type aliases, so any class that defines them directly (e.g. [**IdentityInterpolationBuilder**](classIdentityInterpolationBuilder.md), SplineBuilder1D) satisfies the InterpolationBuilder concept without specialisation.


Specialise this struct to adapt external builders whose alias names differ from the convention (e.g. ddc::SplineBuilder).




**Template parameters:**


* `Builder` The interpolation builder type. 




    
## Public Types Documentation




### typedef basis\_domain\_type 

_The discrete dimension for the interpolation coefficients._ 
```C++
using InterpolationBuilderTraits< Builder >::basis_domain_type =  typename Builder::basis_domain_type;
```




<hr>



### typedef batched\_basis\_idx\_range\_type 

_Batched domain with interpolation\_grid\_type replaced by basis\_domain\_type._ 
```C++
using InterpolationBuilderTraits< Builder >::batched_basis_idx_range_type =  typename Builder::template batched_basis_idx_range_type<IdxRangeBatchedInterpolation>;
```




<hr>



### typedef batched\_derivs\_idx\_range\_type 

_Batched domain with interpolation\_grid\_type replaced by deriv\_type._ 
```C++
using InterpolationBuilderTraits< Builder >::batched_derivs_idx_range_type =  typename Builder::template batched_derivs_idx_range_type<IdxRangeBatchedInterpolation>;
```




<hr>



### typedef data\_type 

_The data type that the data is saved on._ 
```C++
using InterpolationBuilderTraits< Builder >::data_type =  typename Builder::data_type;
```




<hr>



### typedef interpolation\_grid\_type 

_The discrete grid on which interpolation values are given._ 
```C++
using InterpolationBuilderTraits< Builder >::interpolation_grid_type =  typename Builder::interpolation_grid_type;
```




<hr>



### typedef interpolation\_idx\_range\_type 

_The 1D index range for the interpolation mesh._ 
```C++
using InterpolationBuilderTraits< Builder >::interpolation_idx_range_type =  typename Builder::interpolation_idx_range_type;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/interpolation/i_interpolation_builder.hpp`

