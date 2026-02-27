

# Struct InterpolationBuilderTraits&lt; ddc::SplineBuilder&lt; ExecSpace, MemorySpace, BSplines, InterpolationDDim, BcLower, BcUpper, Solver &gt; &gt;

**template &lt;class ExecSpace, class MemorySpace, class BSplines, class InterpolationDDim, ddc::BoundCond BcLower, ddc::BoundCond BcUpper, ddc::SplineSolver Solver&gt;**



[**ClassList**](annotated.md) **>** [**InterpolationBuilderTraits&lt; ddc::SplineBuilder&lt; ExecSpace, MemorySpace, BSplines, InterpolationDDim, BcLower, BcUpper, Solver &gt; &gt;**](structInterpolationBuilderTraits_3_01ddc_1_1SplineBuilder_3_01ExecSpace_00_01MemorySpace_00_01BScfda3229aca4044474c5fb515881d93f.md)



_Specialisation of_ [_**InterpolationBuilderTraits**_](structInterpolationBuilderTraits.md) _for ddc::SplineBuilder._[More...](#detailed-description)

* `#include <i_interpolation_builder.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef typename Builder::bsplines\_type | [**basis\_domain\_type**](#typedef-basis_domain_type)  <br>_The discrete dimension for the B-spline coefficients._  |
| typedef typename Builder::template batched\_spline\_domain\_type&lt; IdxRangeBatchedInterpolation &gt; | [**batched\_basis\_idx\_range\_type**](#typedef-batched_basis_idx_range_type)  <br>_Batched domain with InterpolationDDim replaced by BSplines._  |
| typedef typename Builder::template batched\_derivs\_domain\_type&lt; IdxRangeBatchedInterpolation &gt; | [**batched\_derivs\_idx\_range\_type**](#typedef-batched_derivs_idx_range_type)  <br>_Batched domain with InterpolationDDim replaced by deriv\_type._  |
| typedef double | [**data\_type**](#typedef-data_type)  <br>_The data type that the data is saved on._  |
| typedef typename Builder::interpolation\_discrete\_dimension\_type | [**interpolation\_grid\_type**](#typedef-interpolation_grid_type)  <br>_The discrete grid on which interpolation values are given._  |
| typedef typename Builder::interpolation\_domain\_type | [**interpolation\_idx\_range\_type**](#typedef-interpolation_idx_range_type)  <br>_The 1D index range for the interpolation mesh._  |
















































## Detailed Description


ddc::SplineBuilder uses different alias names from the InterpolationBuilder convention. This specialisation provides the mapping so that ddc::SplineBuilder can be used directly as an InterpolationBuilder without wrapping it in SplineBuilder1D.


Mapping: interpolation\_discrete\_dimension\_type -&gt; interpolation\_grid\_type interpolation\_domain\_type -&gt; interpolation\_idx\_range\_type bsplines\_type -&gt; basis\_domain\_type batched\_spline\_domain\_type&lt;D&gt; -&gt; batched\_basis\_idx\_range\_type&lt;D&gt; batched\_derivs\_domain\_type&lt;D&gt; -&gt; batched\_derivs\_idx\_range\_type&lt;D&gt; 


    
## Public Types Documentation




### typedef basis\_domain\_type 

_The discrete dimension for the B-spline coefficients._ 
```C++
using InterpolationBuilderTraits< ddc::SplineBuilder< ExecSpace, MemorySpace, BSplines, InterpolationDDim, BcLower, BcUpper, Solver > >::basis_domain_type =  typename Builder::bsplines_type;
```




<hr>



### typedef batched\_basis\_idx\_range\_type 

_Batched domain with InterpolationDDim replaced by BSplines._ 
```C++
using InterpolationBuilderTraits< ddc::SplineBuilder< ExecSpace, MemorySpace, BSplines, InterpolationDDim, BcLower, BcUpper, Solver > >::batched_basis_idx_range_type =  typename Builder::template batched_spline_domain_type<IdxRangeBatchedInterpolation>;
```




<hr>



### typedef batched\_derivs\_idx\_range\_type 

_Batched domain with InterpolationDDim replaced by deriv\_type._ 
```C++
using InterpolationBuilderTraits< ddc::SplineBuilder< ExecSpace, MemorySpace, BSplines, InterpolationDDim, BcLower, BcUpper, Solver > >::batched_derivs_idx_range_type =  typename Builder::template batched_derivs_domain_type<IdxRangeBatchedInterpolation>;
```




<hr>



### typedef data\_type 

_The data type that the data is saved on._ 
```C++
using InterpolationBuilderTraits< ddc::SplineBuilder< ExecSpace, MemorySpace, BSplines, InterpolationDDim, BcLower, BcUpper, Solver > >::data_type =  double;
```




<hr>



### typedef interpolation\_grid\_type 

_The discrete grid on which interpolation values are given._ 
```C++
using InterpolationBuilderTraits< ddc::SplineBuilder< ExecSpace, MemorySpace, BSplines, InterpolationDDim, BcLower, BcUpper, Solver > >::interpolation_grid_type =  typename Builder::interpolation_discrete_dimension_type;
```




<hr>



### typedef interpolation\_idx\_range\_type 

_The 1D index range for the interpolation mesh._ 
```C++
using InterpolationBuilderTraits< ddc::SplineBuilder< ExecSpace, MemorySpace, BSplines, InterpolationDDim, BcLower, BcUpper, Solver > >::interpolation_idx_range_type =  typename Builder::interpolation_domain_type;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/interpolation/i_interpolation_builder.hpp`

