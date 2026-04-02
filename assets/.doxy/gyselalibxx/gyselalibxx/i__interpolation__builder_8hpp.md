

# File i\_interpolation\_builder.hpp



[**FileList**](files.md) **>** [**interpolation**](dir_264890e5c091f8c8d7fe1f842870c25e.md) **>** [**i\_interpolation\_builder.hpp**](i__interpolation__builder_8hpp.md)

[Go to the source code of this file](i__interpolation__builder_8hpp_source.md)



* `#include <optional>`
* `#include <type_traits>`
* `#include "ddc_aliases.hpp"`













## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**concepts**](namespaceconcepts.md) <br> |


## Classes

| Type | Name |
| ---: | :--- |
| struct | [**InterpolationBuilderTraits**](structInterpolationBuilderTraits.md) &lt;class Builder&gt;<br>_A traits struct for accessing type aliases of an interpolation builder._  |
| struct | [**InterpolationBuilderTraits&lt; ddc::SplineBuilder&lt; ExecSpace, MemorySpace, BSplines, InterpolationDDim, BcLower, BcUpper, Solver &gt; &gt;**](structInterpolationBuilderTraits_3_01ddc_1_1SplineBuilder_3_01ExecSpace_00_01MemorySpace_00_01BScfda3229aca4044474c5fb515881d93f.md) &lt;class ExecSpace, class MemorySpace, class BSplines, class InterpolationDDim, BcLower, BcUpper, Solver&gt;<br>_Specialisation of_ [_**InterpolationBuilderTraits**_](structInterpolationBuilderTraits.md) _for ddc::SplineBuilder._ |


## Public Types

| Type | Name |
| ---: | :--- |
| typedef ddc::type\_seq\_element\_t&lt; 0, ddc::to\_type\_seq\_t&lt; typename [**InterpolationBuilderTraits**](structInterpolationBuilderTraits.md)&lt; BuilderType &gt;::interpolation\_idx\_range\_type &gt; &gt; | [**interpolation\_grid\_t**](#typedef-interpolation_grid_t)  <br>_The discrete grid on which interpolation values are given (1D builders only)._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
|  auto | [**batched\_basis\_idx\_range**](#function-batched_basis_idx_range) (Builder const & builder, IdxRangeBatchedInterpolation const & batched\_interpolation\_domain) <br>_Get the batched basis index range for a builder._  |




























## Public Types Documentation




### typedef interpolation\_grid\_t 

_The discrete grid on which interpolation values are given (1D builders only)._ 
```C++
using interpolation_grid_t =  ddc::type_seq_element_t< 0, ddc::to_type_seq_t< typename InterpolationBuilderTraits<BuilderType>::interpolation_idx_range_type> >;
```




<hr>
## Public Functions Documentation




### function batched\_basis\_idx\_range 

_Get the batched basis index range for a builder._ 
```C++
template<class Builder, class IdxRangeBatchedInterpolation>
auto batched_basis_idx_range (
    Builder const & builder,
    IdxRangeBatchedInterpolation const & batched_interpolation_domain
) 
```



Dispatches to batched\_basis\_idx\_range if available (e.g. [**IdentityInterpolationBuilder**](classIdentityInterpolationBuilder.md)), otherwise falls back to batched\_spline\_domain (e.g. ddc::SplineBuilder).




**Parameters:**


* `builder` The interpolation builder. 
* `batched_interpolation_domain` The batched interpolation domain. 



**Returns:**

The batched basis index range. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/interpolation/i_interpolation_builder.hpp`

