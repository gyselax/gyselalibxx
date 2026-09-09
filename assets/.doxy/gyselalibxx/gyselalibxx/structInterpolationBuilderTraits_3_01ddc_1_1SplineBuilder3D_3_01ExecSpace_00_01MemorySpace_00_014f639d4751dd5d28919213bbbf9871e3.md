

# Struct InterpolationBuilderTraits&lt; ddc::SplineBuilder3D&lt; ExecSpace, MemorySpace, BSpline1, BSpline2, BSpline3, InterpolationDDim1, InterpolationDDim2, InterpolationDDim3, BcLower1, BcUpper1, BcLower2, BcUpper2, BcLower3, BcUpper3, Solver &gt; &gt;

**template &lt;class ExecSpace, class MemorySpace, class BSpline1, class BSpline2, class BSpline3, class InterpolationDDim1, class InterpolationDDim2, class InterpolationDDim3, ddc::SplineBuilderClosure BcLower1, ddc::SplineBuilderClosure BcUpper1, ddc::SplineBuilderClosure BcLower2, ddc::SplineBuilderClosure BcUpper2, ddc::SplineBuilderClosure BcLower3, ddc::SplineBuilderClosure BcUpper3, ddc::SplineSolver Solver&gt;**



[**ClassList**](annotated.md) **>** [**InterpolationBuilderTraits&lt; ddc::SplineBuilder3D&lt; ExecSpace, MemorySpace, BSpline1, BSpline2, BSpline3, InterpolationDDim1, InterpolationDDim2, InterpolationDDim3, BcLower1, BcUpper1, BcLower2, BcUpper2, BcLower3, BcUpper3, Solver &gt; &gt;**](structInterpolationBuilderTraits_3_01ddc_1_1SplineBuilder3D_3_01ExecSpace_00_01MemorySpace_00_014f639d4751dd5d28919213bbbf9871e3.md)



_Specialisation of_ [_**InterpolationBuilderTraits**_](structInterpolationBuilderTraits.md) _for ddc::SplineBuilder3D._[More...](#detailed-description)

* `#include <i_interpolation_builder.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef typename Builder::template batched\_spline\_domain\_type&lt; IdxRangeBatchedInterpolation &gt; | [**batched\_basis\_idx\_range\_type**](#typedef-batched_basis_idx_range_type)  <br>_Batched domain with InterpolationDDim replaced by BSplines._  |
| typedef typename Builder::template batched\_derivs\_domain\_type&lt; IdxRangeBatchedInterpolation &gt; | [**batched\_derivs\_idx\_range\_type**](#typedef-batched_derivs_idx_range_type)  <br>_Batched domain with InterpolationDDim replaced by deriv\_type._  |
| typedef IdxRange&lt; typename Builder::bsplines\_type1, typename Builder::bsplines\_type2, typename Builder::bsplines\_type3 &gt; | [**coeff\_idx\_range\_type**](#typedef-coeff_idx_range_type)  <br>_The discrete dimension for the B-spline coefficients._  |
| typedef double | [**data\_type**](#typedef-data_type)  <br>_The data type that the data is saved on._  |
| typedef typename Builder::interpolation\_domain\_type | [**interpolation\_idx\_range\_type**](#typedef-interpolation_idx_range_type)  <br>_The discrete grid on which interpolation values are given._  |






















## Public Static Functions

| Type | Name |
| ---: | :--- |
|  constexpr std::size\_t | [**rank**](#function-rank) () <br>_The number of interpolation dimensions (always 1 for SplineBuilder)._  |


























## Detailed Description


ddc::SplineBuilder3D uses different alias names from the InterpolationBuilder convention. This specialisation provides the mapping so that ddc::SplineBuilder3D can be used directly as an InterpolationBuilder without wrapping it.


Mapping: interpolation\_discrete\_dimension\_type -&gt; interpolation\_grid\_type interpolation\_domain\_type -&gt; interpolation\_idx\_range\_type bsplines\_type -&gt; basis\_domain\_type batched\_spline\_domain\_type&lt;D&gt; -&gt; batched\_basis\_idx\_range\_type&lt;D&gt; batched\_derivs\_domain\_type&lt;D&gt; -&gt; batched\_derivs\_idx\_range\_type&lt;D&gt; 


    
## Public Types Documentation




### typedef batched\_basis\_idx\_range\_type 

_Batched domain with InterpolationDDim replaced by BSplines._ 
```C++
using InterpolationBuilderTraits< ddc::SplineBuilder3D< ExecSpace, MemorySpace, BSpline1, BSpline2, BSpline3, InterpolationDDim1, InterpolationDDim2, InterpolationDDim3, BcLower1, BcUpper1, BcLower2, BcUpper2, BcLower3, BcUpper3, Solver > >::batched_basis_idx_range_type =  typename Builder::template batched_spline_domain_type<IdxRangeBatchedInterpolation>;
```




<hr>



### typedef batched\_derivs\_idx\_range\_type 

_Batched domain with InterpolationDDim replaced by deriv\_type._ 
```C++
using InterpolationBuilderTraits< ddc::SplineBuilder3D< ExecSpace, MemorySpace, BSpline1, BSpline2, BSpline3, InterpolationDDim1, InterpolationDDim2, InterpolationDDim3, BcLower1, BcUpper1, BcLower2, BcUpper2, BcLower3, BcUpper3, Solver > >::batched_derivs_idx_range_type =  typename Builder::template batched_derivs_domain_type<IdxRangeBatchedInterpolation>;
```




<hr>



### typedef coeff\_idx\_range\_type 

_The discrete dimension for the B-spline coefficients._ 
```C++
using InterpolationBuilderTraits< ddc::SplineBuilder3D< ExecSpace, MemorySpace, BSpline1, BSpline2, BSpline3, InterpolationDDim1, InterpolationDDim2, InterpolationDDim3, BcLower1, BcUpper1, BcLower2, BcUpper2, BcLower3, BcUpper3, Solver > >::coeff_idx_range_type =  IdxRange< typename Builder::bsplines_type1, typename Builder::bsplines_type2, typename Builder::bsplines_type3>;
```




<hr>



### typedef data\_type 

_The data type that the data is saved on._ 
```C++
using InterpolationBuilderTraits< ddc::SplineBuilder3D< ExecSpace, MemorySpace, BSpline1, BSpline2, BSpline3, InterpolationDDim1, InterpolationDDim2, InterpolationDDim3, BcLower1, BcUpper1, BcLower2, BcUpper2, BcLower3, BcUpper3, Solver > >::data_type =  double;
```




<hr>



### typedef interpolation\_idx\_range\_type 

_The discrete grid on which interpolation values are given._ 
```C++
using InterpolationBuilderTraits< ddc::SplineBuilder3D< ExecSpace, MemorySpace, BSpline1, BSpline2, BSpline3, InterpolationDDim1, InterpolationDDim2, InterpolationDDim3, BcLower1, BcUpper1, BcLower2, BcUpper2, BcLower3, BcUpper3, Solver > >::interpolation_idx_range_type =  typename Builder::interpolation_domain_type;
```



The 1D index range for the interpolation mesh. 


        

<hr>
## Public Static Functions Documentation




### function rank 

_The number of interpolation dimensions (always 1 for SplineBuilder)._ 
```C++
static inline constexpr std::size_t InterpolationBuilderTraits< ddc::SplineBuilder3D< ExecSpace, MemorySpace, BSpline1, BSpline2, BSpline3, InterpolationDDim1, InterpolationDDim2, InterpolationDDim3, BcLower1, BcUpper1, BcLower2, BcUpper2, BcLower3, BcUpper3, Solver > >::rank () 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/interpolation/i_interpolation_builder.hpp`

