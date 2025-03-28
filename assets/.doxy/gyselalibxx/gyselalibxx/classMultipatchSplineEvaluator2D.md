

# Class MultipatchSplineEvaluator2D

**template &lt;class ExecSpace, class MemorySpace, template&lt; typename P &gt; typename BSpline1OnPatch, template&lt; typename P &gt; typename BSpline2OnPatch, template&lt; typename P &gt; typename Grid1OnPatch, template&lt; typename P &gt; typename Grid2OnPatch, class ExtrapolationRule, template&lt; typename P &gt; typename ValuesOnPatch, class PatchLocator, class... Patches&gt;**



[**ClassList**](annotated.md) **>** [**MultipatchSplineEvaluator2D**](classMultipatchSplineEvaluator2D.md)



_A class to evaluate all the splines of all the patches at once._ [More...](#detailed-description)

* `#include <multipatch_spline_evaluator_2d.hpp>`















## Classes

| Type | Name |
| ---: | :--- |
| struct | [**eval\_deriv\_type**](structMultipatchSplineEvaluator2D_1_1eval__deriv__type.md) <br>_Tag to indicate that derivative of the spline should be evaluated._  |
| struct | [**eval\_type**](structMultipatchSplineEvaluator2D_1_1eval__type.md) <br>_Tag to indicate that the value of the spline should be evaluated._  |


## Public Types

| Type | Name |
| ---: | :--- |
| typedef DConstField&lt; [**spline\_idx\_range\_type**](classMultipatchSplineEvaluator2D.md#typedef-spline_idx_range_type)&lt; [**Patch**](structPatch.md) &gt;, MemorySpace &gt; | [**SplineCoeffOnPatch**](#typedef-splinecoeffonpatch)  <br>_Type for_ [_**MultipatchType**_](classMultipatchType.md) _: A field of 2D spline coefficients for a non-batched spline defined on both of the_[_**Patch**_](structPatch.md) _'s logical dimensions. Needed public for functions on GPU._ |
| typedef typename ValuesOnPatch&lt; [**Patch**](structPatch.md) &gt;::discrete\_domain\_type | [**batched\_evaluation\_idx\_range\_type**](#typedef-batched_evaluation_idx_range_type)  <br>_The type of the whole domain representing evaluation points._  |
| typedef BSpline1OnPatch&lt; [**Patch**](structPatch.md) &gt; | [**bsplines\_type1**](#typedef-bsplines_type1)  <br>_The discrete dimension representing the B-splines along first dimension._  |
| typedef BSpline2OnPatch&lt; [**Patch**](structPatch.md) &gt; | [**bsplines\_type2**](#typedef-bsplines_type2)  <br>_The discrete dimension representing the B-splines along second dimension._  |
| typedef typename BSpline1OnPatch&lt; [**Patch**](structPatch.md) &gt;::continuous\_dimension\_type | [**continuous\_dimension\_type1**](#typedef-continuous_dimension_type1)  <br>_The type of the first evaluation continuous dimension used by this class._  |
| typedef typename BSpline2OnPatch&lt; [**Patch**](structPatch.md) &gt;::continuous\_dimension\_type | [**continuous\_dimension\_type2**](#typedef-continuous_dimension_type2)  <br>_The type of the second evaluation continuous dimension used by this class._  |
| typedef Grid1OnPatch&lt; [**Patch**](structPatch.md) &gt; | [**evaluation\_discrete\_dimension\_type1**](#typedef-evaluation_discrete_dimension_type1)  <br>_The type of the first discrete dimension of interest used by this class._  |
| typedef Grid2OnPatch&lt; [**Patch**](structPatch.md) &gt; | [**evaluation\_discrete\_dimension\_type2**](#typedef-evaluation_discrete_dimension_type2)  <br>_The type of the second discrete dimension of interest used by this class._  |
| typedef IdxRange&lt; [**evaluation\_discrete\_dimension\_type1**](classMultipatchSplineEvaluator2D.md#typedef-evaluation_discrete_dimension_type1)&lt; [**Patch**](structPatch.md) &gt;, [**evaluation\_discrete\_dimension\_type2**](classMultipatchSplineEvaluator2D.md#typedef-evaluation_discrete_dimension_type2)&lt; [**Patch**](structPatch.md) &gt; &gt; | [**evaluation\_idx\_range\_type**](#typedef-evaluation_idx_range_type)  <br>_The type of the domain for the 2D evaluation mesh used by this class._  |
| typedef IdxRange&lt; [**evaluation\_discrete\_dimension\_type1**](classMultipatchSplineEvaluator2D.md#typedef-evaluation_discrete_dimension_type1)&lt; [**Patch**](structPatch.md) &gt; &gt; | [**evaluation\_idx\_range\_type1**](#typedef-evaluation_idx_range_type1)  <br>_The type of the domain for the 1D evaluation mesh along first dimension used by this class._  |
| typedef IdxRange&lt; [**evaluation\_discrete\_dimension\_type2**](classMultipatchSplineEvaluator2D.md#typedef-evaluation_discrete_dimension_type2)&lt; [**Patch**](structPatch.md) &gt; &gt; | [**evaluation\_idx\_range\_type2**](#typedef-evaluation_idx_range_type2)  <br>_The type of the domain for the 1D evaluation mesh along second dimension used by this class._  |
| typedef ExecSpace | [**exec\_space**](#typedef-exec_space)  <br>_The type of the Kokkos execution space used by this class._  |
| typedef MemorySpace | [**memory\_space**](#typedef-memory_space)  <br>_The type of the Kokkos memory space used by this class._  |
| typedef IdxRange&lt; [**bsplines\_type1**](classMultipatchSplineEvaluator2D.md#typedef-bsplines_type1)&lt; [**Patch**](structPatch.md) &gt;, [**bsplines\_type2**](classMultipatchSplineEvaluator2D.md#typedef-bsplines_type2)&lt; [**Patch**](structPatch.md) &gt; &gt; | [**spline\_idx\_range\_type**](#typedef-spline_idx_range_type)  <br>_The type of the 2D spline domain corresponding to the dimensions of interest._  |
| typedef IdxRange&lt; [**bsplines\_type1**](classMultipatchSplineEvaluator2D.md#typedef-bsplines_type1)&lt; [**Patch**](structPatch.md) &gt; &gt; | [**spline\_idx\_range\_type1**](#typedef-spline_idx_range_type1)  <br>_The type of the 1D spline domain corresponding to the first dimension of interest._  |
| typedef IdxRange&lt; [**bsplines\_type2**](classMultipatchSplineEvaluator2D.md#typedef-bsplines_type2)&lt; [**Patch**](structPatch.md) &gt; &gt; | [**spline\_idx\_range\_type2**](#typedef-spline_idx_range_type2)  <br>_The type of the 1D spline domain corresponding to the second dimension of interest._  |






## Public Static Attributes

| Type | Name |
| ---: | :--- |
|  constexpr std::size\_t | [**n\_patches**](#variable-n_patches)   = `ddc::type\_seq\_size\_v&lt;PatchOrdering&gt;`<br>_The number of patches._  |














## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**MultipatchSplineEvaluator2D**](#function-multipatchsplineevaluator2d) (PatchLocator const & patch\_locator, ExtrapolationRule const & extrapolation\_rule) <br>_Instantiate a_ [_**MultipatchSplineEvaluator2D**_](classMultipatchSplineEvaluator2D.md) _._ |
|  void | [**apply\_evaluator**](#function-apply_evaluator) (ValuesOnPatch&lt; StoringPatch &gt; const & patch\_values, CoordConstFieldOnPatch&lt; StoringPatch &gt; const & patch\_coords, [**MultipatchSplineCoeff**](classMultipatchField.md) const patches\_splines) const<br>_Compute the values or the derivatives of a given patch at the coordinates defined on the given patch. Needed public for functions on GPU._  |
|  void | [**apply\_integrate**](#function-apply_integrate) (double & integral, [**SplineCoeffOnPatch**](classMultipatchSplineEvaluator2D.md#typedef-splinecoeffonpatch)&lt; [**Patch**](structPatch.md) &gt; const & spline\_coef) const<br>_Integrate the spline defined on the given patch._  |
|  KOKKOS\_FUNCTION double | [**deriv**](#function-deriv) (Coord&lt; Dim1, Dim2 &gt; const & coord\_eval, [**MultipatchSplineCoeff**](classMultipatchField.md) const & patches\_splines) const<br>_Differentiate 2D splines (described by their spline coefficients) at a given coordinate along a specified dimension of interest._  |
|  KOKKOS\_FUNCTION double | [**deriv\_1\_and\_2**](#function-deriv_1_and_2-12) (Coord const & coord\_eval, [**MultipatchSplineCoeff**](classMultipatchField.md) const & patches\_splines) const<br>_Cross-differentiate 2D splines (described by their spline coefficients) at a given coordinate._  |
|  void | [**deriv\_1\_and\_2**](#function-deriv_1_and_2-22) ([**MultipatchValues**](classMultipatchField.md) const & patches\_deriv\_12, [**MultipatchCoordField**](classMultipatchField.md) const & patches\_coords, [**MultipatchSplineCoeff**](classMultipatchField.md) const & patches\_splines) const<br>_Cross-differentiate 2D splines (described by their spline coefficients) on a meshes._  |
|  KOKKOS\_FUNCTION double | [**deriv\_dim\_1**](#function-deriv_dim_1-12) (Coord const & coord\_eval, [**MultipatchSplineCoeff**](classMultipatchField.md) const & patches\_splines) const<br>_Differentiate 2D splines (described by their spline coefficients) at a given coordinate along first dimension of interest._  |
|  void | [**deriv\_dim\_1**](#function-deriv_dim_1-22) ([**MultipatchValues**](classMultipatchField.md) const & patches\_deriv\_1, [**MultipatchCoordField**](classMultipatchField.md) const & patches\_coords, [**MultipatchSplineCoeff**](classMultipatchField.md) const & patches\_splines) const<br>_Differentiate 2D splines (described by their spline coefficients) on a meshes along first dimension of interest._  |
|  KOKKOS\_FUNCTION double | [**deriv\_dim\_2**](#function-deriv_dim_2-12) (Coord const & coord\_eval, [**MultipatchSplineCoeff**](classMultipatchField.md) const & patches\_splines) const<br>_Differentiate 2D splines (described by their spline coefficients) at a given coordinate along second dimension of interest._  |
|  void | [**deriv\_dim\_2**](#function-deriv_dim_2-22) ([**MultipatchValues**](classMultipatchField.md) const & patches\_deriv\_2, [**MultipatchCoordField**](classMultipatchField.md) const & patches\_coords, [**MultipatchSplineCoeff**](classMultipatchField.md) const & patches\_splines) const<br>_Differentiate 2D splines (described by their spline coefficients) on a meshes along second dimension of interest._  |
|  void | [**integrate**](#function-integrate) (DKokkosView\_h&lt; [**n\_patches**](classMultipatchSplineEvaluator2D.md#variable-n_patches) &gt; const & integrals, [**MultipatchSplineCoeff**](classMultipatchField.md) const & patches\_splines) const<br>_Integration of splines (described by their spline coefficients)._  |
|  KOKKOS\_FUNCTION double | [**operator()**](#function-operator) (Coord const coord\_eval, [**MultipatchSplineCoeff**](classMultipatchField.md) const & patches\_splines) const<br>_Evaluate the 2D splines (described by their spline coefficients) at a given coordinate._  |
|  void | [**operator()**](#function-operator_1) ([**MultipatchValues**](classMultipatchField.md) const & patches\_values, [**MultipatchCoordField**](classMultipatchField.md) const & patches\_coords, [**MultipatchSplineCoeff**](classMultipatchField.md) const & patches\_splines) const<br>_Evaluate 2D splines (described by their spline coefficients) on meshes._  |
|   | [**~MultipatchSplineEvaluator2D**](#function-multipatchsplineevaluator2d) () = default<br> |




























## Detailed Description


This class computes the evaluations of the splines defined on each patch at given coordinates. The class does not need to instantiate in advance individual spline evaluators. The coordinates are stored in fields on each patch. On a given storing patch, we do not enforce that a coordinate is physically located on this storing patch. These fields of coordinates can be the result of characteristics equations solving.


This function is useful to avoid calling all the spline evaluators individually, especially in a multipatch geometry with several patches. It also manages to evaluate the right spline at a coordinate physically located on another patch than the storing patch.


Additionally, methods to compute the first derivatives and cross-derivatives are implemented.




**Warning:**

This operator does not work on batched domain.




**Template parameters:**


* `ExecSpace` The space (CPU/GPU) where the calculations are carried out. 
* `MemorySpace` The space (CPU/GPU) where the coefficients and values are stored. 
* `BSpline1OnPatch` A type alias which provides the first BSpline type along which the splines are built template on the [**Patch**](structPatch.md). 
* `BSpline2OnPatch` A type alias which provides the second BSpline type along which the splines are built template on the [**Patch**](structPatch.md). 
* `Grid1OnPatch` A type alias which provides the first Grid type along which the interpolation points of the splines are found template on the [**Patch**](structPatch.md). 
* `Grid2OnPatch` A type alias which provides the second Grid type along which the interpolation points of the splines are found template on the [**Patch**](structPatch.md). 
* `ExtrapolationRule` The extrapolation rule type for outside of the global domain. 
* `ValuesOnPatch` A Field type storing the evaluated values of the splines. Template on the [**Patch**](structPatch.md). 
* `PatchLocator` A operator that finds the patch where a given coordinate is physically located. 
* `Patches` A variadic template of all the patches. 




    
## Public Types Documentation




### typedef SplineCoeffOnPatch 

_Type for_ [_**MultipatchType**_](classMultipatchType.md) _: A field of 2D spline coefficients for a non-batched spline defined on both of the_[_**Patch**_](structPatch.md) _'s logical dimensions. Needed public for functions on GPU._
```C++
using MultipatchSplineEvaluator2D< ExecSpace, MemorySpace, BSpline1OnPatch, BSpline2OnPatch, Grid1OnPatch, Grid2OnPatch, ExtrapolationRule, ValuesOnPatch, PatchLocator, Patches >::SplineCoeffOnPatch =  DConstField<spline_idx_range_type<Patch>, MemorySpace>;
```




<hr>



### typedef batched\_evaluation\_idx\_range\_type 

_The type of the whole domain representing evaluation points._ 
```C++
using MultipatchSplineEvaluator2D< ExecSpace, MemorySpace, BSpline1OnPatch, BSpline2OnPatch, Grid1OnPatch, Grid2OnPatch, ExtrapolationRule, ValuesOnPatch, PatchLocator, Patches >::batched_evaluation_idx_range_type =  typename ValuesOnPatch<Patch>::discrete_domain_type;
```





**Template parameters:**


* [**Patch**](structPatch.md) [**Patch**](structPatch.md) type. 




        

<hr>



### typedef bsplines\_type1 

_The discrete dimension representing the B-splines along first dimension._ 
```C++
using MultipatchSplineEvaluator2D< ExecSpace, MemorySpace, BSpline1OnPatch, BSpline2OnPatch, Grid1OnPatch, Grid2OnPatch, ExtrapolationRule, ValuesOnPatch, PatchLocator, Patches >::bsplines_type1 =  BSpline1OnPatch<Patch>;
```





**Template parameters:**


* [**Patch**](structPatch.md) [**Patch**](structPatch.md) type. 




        

<hr>



### typedef bsplines\_type2 

_The discrete dimension representing the B-splines along second dimension._ 
```C++
using MultipatchSplineEvaluator2D< ExecSpace, MemorySpace, BSpline1OnPatch, BSpline2OnPatch, Grid1OnPatch, Grid2OnPatch, ExtrapolationRule, ValuesOnPatch, PatchLocator, Patches >::bsplines_type2 =  BSpline2OnPatch<Patch>;
```





**Template parameters:**


* [**Patch**](structPatch.md) [**Patch**](structPatch.md) type. 




        

<hr>



### typedef continuous\_dimension\_type1 

_The type of the first evaluation continuous dimension used by this class._ 
```C++
using MultipatchSplineEvaluator2D< ExecSpace, MemorySpace, BSpline1OnPatch, BSpline2OnPatch, Grid1OnPatch, Grid2OnPatch, ExtrapolationRule, ValuesOnPatch, PatchLocator, Patches >::continuous_dimension_type1 =  typename BSpline1OnPatch<Patch>::continuous_dimension_type;
```





**Template parameters:**


* [**Patch**](structPatch.md) [**Patch**](structPatch.md) type. 




        

<hr>



### typedef continuous\_dimension\_type2 

_The type of the second evaluation continuous dimension used by this class._ 
```C++
using MultipatchSplineEvaluator2D< ExecSpace, MemorySpace, BSpline1OnPatch, BSpline2OnPatch, Grid1OnPatch, Grid2OnPatch, ExtrapolationRule, ValuesOnPatch, PatchLocator, Patches >::continuous_dimension_type2 =  typename BSpline2OnPatch<Patch>::continuous_dimension_type;
```





**Template parameters:**


* [**Patch**](structPatch.md) [**Patch**](structPatch.md) type. 




        

<hr>



### typedef evaluation\_discrete\_dimension\_type1 

_The type of the first discrete dimension of interest used by this class._ 
```C++
using MultipatchSplineEvaluator2D< ExecSpace, MemorySpace, BSpline1OnPatch, BSpline2OnPatch, Grid1OnPatch, Grid2OnPatch, ExtrapolationRule, ValuesOnPatch, PatchLocator, Patches >::evaluation_discrete_dimension_type1 =  Grid1OnPatch<Patch>;
```





**Template parameters:**


* [**Patch**](structPatch.md) [**Patch**](structPatch.md) type. 




        

<hr>



### typedef evaluation\_discrete\_dimension\_type2 

_The type of the second discrete dimension of interest used by this class._ 
```C++
using MultipatchSplineEvaluator2D< ExecSpace, MemorySpace, BSpline1OnPatch, BSpline2OnPatch, Grid1OnPatch, Grid2OnPatch, ExtrapolationRule, ValuesOnPatch, PatchLocator, Patches >::evaluation_discrete_dimension_type2 =  Grid2OnPatch<Patch>;
```





**Template parameters:**


* [**Patch**](structPatch.md) [**Patch**](structPatch.md) type. 




        

<hr>



### typedef evaluation\_idx\_range\_type 

_The type of the domain for the 2D evaluation mesh used by this class._ 
```C++
using MultipatchSplineEvaluator2D< ExecSpace, MemorySpace, BSpline1OnPatch, BSpline2OnPatch, Grid1OnPatch, Grid2OnPatch, ExtrapolationRule, ValuesOnPatch, PatchLocator, Patches >::evaluation_idx_range_type =  IdxRange< evaluation_discrete_dimension_type1<Patch>, evaluation_discrete_dimension_type2<Patch> >;
```





**Template parameters:**


* [**Patch**](structPatch.md) [**Patch**](structPatch.md) type. 




        

<hr>



### typedef evaluation\_idx\_range\_type1 

_The type of the domain for the 1D evaluation mesh along first dimension used by this class._ 
```C++
using MultipatchSplineEvaluator2D< ExecSpace, MemorySpace, BSpline1OnPatch, BSpline2OnPatch, Grid1OnPatch, Grid2OnPatch, ExtrapolationRule, ValuesOnPatch, PatchLocator, Patches >::evaluation_idx_range_type1 =  IdxRange<evaluation_discrete_dimension_type1<Patch> >;
```





**Template parameters:**


* [**Patch**](structPatch.md) [**Patch**](structPatch.md) type. 




        

<hr>



### typedef evaluation\_idx\_range\_type2 

_The type of the domain for the 1D evaluation mesh along second dimension used by this class._ 
```C++
using MultipatchSplineEvaluator2D< ExecSpace, MemorySpace, BSpline1OnPatch, BSpline2OnPatch, Grid1OnPatch, Grid2OnPatch, ExtrapolationRule, ValuesOnPatch, PatchLocator, Patches >::evaluation_idx_range_type2 =  IdxRange<evaluation_discrete_dimension_type2<Patch> >;
```





**Template parameters:**


* [**Patch**](structPatch.md) [**Patch**](structPatch.md) type. 




        

<hr>



### typedef exec\_space 

_The type of the Kokkos execution space used by this class._ 
```C++
using MultipatchSplineEvaluator2D< ExecSpace, MemorySpace, BSpline1OnPatch, BSpline2OnPatch, Grid1OnPatch, Grid2OnPatch, ExtrapolationRule, ValuesOnPatch, PatchLocator, Patches >::exec_space =  ExecSpace;
```




<hr>



### typedef memory\_space 

_The type of the Kokkos memory space used by this class._ 
```C++
using MultipatchSplineEvaluator2D< ExecSpace, MemorySpace, BSpline1OnPatch, BSpline2OnPatch, Grid1OnPatch, Grid2OnPatch, ExtrapolationRule, ValuesOnPatch, PatchLocator, Patches >::memory_space =  MemorySpace;
```




<hr>



### typedef spline\_idx\_range\_type 

_The type of the 2D spline domain corresponding to the dimensions of interest._ 
```C++
using MultipatchSplineEvaluator2D< ExecSpace, MemorySpace, BSpline1OnPatch, BSpline2OnPatch, Grid1OnPatch, Grid2OnPatch, ExtrapolationRule, ValuesOnPatch, PatchLocator, Patches >::spline_idx_range_type =  IdxRange<bsplines_type1<Patch>, bsplines_type2<Patch> >;
```





**Template parameters:**


* [**Patch**](structPatch.md) [**Patch**](structPatch.md) type. 




        

<hr>



### typedef spline\_idx\_range\_type1 

_The type of the 1D spline domain corresponding to the first dimension of interest._ 
```C++
using MultipatchSplineEvaluator2D< ExecSpace, MemorySpace, BSpline1OnPatch, BSpline2OnPatch, Grid1OnPatch, Grid2OnPatch, ExtrapolationRule, ValuesOnPatch, PatchLocator, Patches >::spline_idx_range_type1 =  IdxRange<bsplines_type1<Patch> >;
```





**Template parameters:**


* [**Patch**](structPatch.md) [**Patch**](structPatch.md) type. 




        

<hr>



### typedef spline\_idx\_range\_type2 

_The type of the 1D spline domain corresponding to the second dimension of interest._ 
```C++
using MultipatchSplineEvaluator2D< ExecSpace, MemorySpace, BSpline1OnPatch, BSpline2OnPatch, Grid1OnPatch, Grid2OnPatch, ExtrapolationRule, ValuesOnPatch, PatchLocator, Patches >::spline_idx_range_type2 =  IdxRange<bsplines_type2<Patch> >;
```





**Template parameters:**


* [**Patch**](structPatch.md) [**Patch**](structPatch.md) type. 




        

<hr>
## Public Static Attributes Documentation




### variable n\_patches 

_The number of patches._ 
```C++
constexpr std::size_t MultipatchSplineEvaluator2D< ExecSpace, MemorySpace, BSpline1OnPatch, BSpline2OnPatch, Grid1OnPatch, Grid2OnPatch, ExtrapolationRule, ValuesOnPatch, PatchLocator, Patches >::n_patches;
```




<hr>
## Public Functions Documentation




### function MultipatchSplineEvaluator2D 

_Instantiate a_ [_**MultipatchSplineEvaluator2D**_](classMultipatchSplineEvaluator2D.md) _._
```C++
inline MultipatchSplineEvaluator2D::MultipatchSplineEvaluator2D (
    PatchLocator const & patch_locator,
    ExtrapolationRule const & extrapolation_rule
) 
```





**Parameters:**


* `patch_locator` An operator to locate a coordinate. The mapping stored in this class has to be invertible. 
* `extrapolation_rule` The extrapolation rules. 




        

<hr>



### function apply\_evaluator 

_Compute the values or the derivatives of a given patch at the coordinates defined on the given patch. Needed public for functions on GPU._ 
```C++
template<class EvalType1, class EvalType2, class StoringPatch>
inline void MultipatchSplineEvaluator2D::apply_evaluator (
    ValuesOnPatch< StoringPatch > const & patch_values,
    CoordConstFieldOnPatch< StoringPatch > const & patch_coords,
    MultipatchSplineCoeff const patches_splines
) const
```





**Template parameters:**


* `EvalType1` Evaluation type: either [**eval\_type**](structMultipatchSplineEvaluator2D_1_1eval__type.md) or [**eval\_deriv\_type**](structMultipatchSplineEvaluator2D_1_1eval__deriv__type.md). 
* `EvalType2` Evaluation type: either [**eval\_type**](structMultipatchSplineEvaluator2D_1_1eval__type.md) or [**eval\_deriv\_type**](structMultipatchSplineEvaluator2D_1_1eval__deriv__type.md). 
* `StoringPatch` [**Patch**](structPatch.md) type where the given coordinates are stored. They are not especially physically located on this patch. 



**Parameters:**


* `patch_values` Field of values of the function or derivative. 
* `patch_coords` ConstField of coordinates defined on the StoringPatch and where we want to evaluate the function or derivative. 
* `patches_splines` [**MultipatchType**](classMultipatchType.md) of spline coefficients of the splines on every patches. 




        

<hr>



### function apply\_integrate 

_Integrate the spline defined on the given patch._ 
```C++
template<class Patch>
inline void MultipatchSplineEvaluator2D::apply_integrate (
    double & integral,
    SplineCoeffOnPatch < Patch > const & spline_coef
) const
```





**Template parameters:**


* [**Patch**](structPatch.md) [**Patch**](structPatch.md) type where the integration of the spline is computed. 



**Parameters:**


* `integral` Double, value of the integral of the spline on the given [**Patch**](structPatch.md). 
* `spline_coef` ConstField of spline coefficients of the spline on the given [**Patch**](structPatch.md). 




        

<hr>



### function deriv 

_Differentiate 2D splines (described by their spline coefficients) at a given coordinate along a specified dimension of interest._ 
```C++
template<class InterestDim, class Dim1, class Dim2>
inline KOKKOS_FUNCTION double MultipatchSplineEvaluator2D::deriv (
    Coord< Dim1, Dim2 > const & coord_eval,
    MultipatchSplineCoeff const & patches_splines
) const
```



See MultipatchSplineEvaluatorOperator.




**Warning:**

The derivative cannot be computed outside of the domain.




**Template parameters:**


* `InterestDim` Dimension of StoringPatch along which differentiation is performed. 
* `StoringPatch` [**Patch**](structPatch.md) type where the given coordinate is stored. It does not mean that the coordinate is physically located on the patch. 



**Parameters:**


* `coord_eval` The coordinate where the spline is differentiated. 
* `patches_splines` A [**MultipatchType**](classMultipatchType.md) of DField storing the 2D spline coefficients.



**Returns:**

The derivative of the spline at the desired coordinate on the right patch. 





        

<hr>



### function deriv\_1\_and\_2 [1/2]

_Cross-differentiate 2D splines (described by their spline coefficients) at a given coordinate._ 
```C++
template<class Coord>
inline KOKKOS_FUNCTION double MultipatchSplineEvaluator2D::deriv_1_and_2 (
    Coord const & coord_eval,
    MultipatchSplineCoeff const & patches_splines
) const
```



See MultipatchSplineEvaluatorOperator.




**Warning:**

The derivative cannot be computed outside of the domain.




**Template parameters:**


* `StoringPatch` [**Patch**](structPatch.md) type where the given coordinate is stored. It does not mean that the coordinate is physically located on the patch. 



**Parameters:**


* `coord_eval` The coordinate where the spline is differentiated. 
* `patches_splines` A [**MultipatchType**](classMultipatchType.md) of DField storing the 2D spline coefficients.



**Returns:**

The derivative of the spline at the desired coordinate on the right patch. 





        

<hr>



### function deriv\_1\_and\_2 [2/2]

_Cross-differentiate 2D splines (described by their spline coefficients) on a meshes._ 
```C++
inline void MultipatchSplineEvaluator2D::deriv_1_and_2 (
    MultipatchValues const & patches_deriv_12,
    MultipatchCoordField const & patches_coords,
    MultipatchSplineCoeff const & patches_splines
) const
```



See MultipatchSplineEvaluatorOperator.




**Warning:**

The derivatives cannot be computed outside of the domain.




**Parameters:**


* `patches_deriv_12` A [**MultipatchType**](classMultipatchType.md) of DField to store the cross-derivatives of the splines at the given coordinates. 
* `patches_coords` A [**MultipatchType**](classMultipatchType.md) of Field of Coordinate storing the coordinates of the meshes. 
* `patches_splines` A [**MultipatchType**](classMultipatchType.md) of DField storing the 2D spline coefficients. 




        

<hr>



### function deriv\_dim\_1 [1/2]

_Differentiate 2D splines (described by their spline coefficients) at a given coordinate along first dimension of interest._ 
```C++
template<class Coord>
inline KOKKOS_FUNCTION double MultipatchSplineEvaluator2D::deriv_dim_1 (
    Coord const & coord_eval,
    MultipatchSplineCoeff const & patches_splines
) const
```



See MultipatchSplineEvaluatorOperator.




**Warning:**

The derivative cannot be computed outside of the domain. The coordinate do not still have to be defined on the right patch.




**Template parameters:**


* `StoringPatch` [**Patch**](structPatch.md) type where the given coordinate is stored. It does not mean that the coordinate is physically located on the patch. 



**Parameters:**


* `coord_eval` The coordinate where the spline is differentiated. 
* `patches_splines` A [**MultipatchType**](classMultipatchType.md) of DField storing the 2D spline coefficients.



**Returns:**

The derivative of the spline at the desired coordinate on the right patch. 





        

<hr>



### function deriv\_dim\_1 [2/2]

_Differentiate 2D splines (described by their spline coefficients) on a meshes along first dimension of interest._ 
```C++
inline void MultipatchSplineEvaluator2D::deriv_dim_1 (
    MultipatchValues const & patches_deriv_1,
    MultipatchCoordField const & patches_coords,
    MultipatchSplineCoeff const & patches_splines
) const
```



See MultipatchSplineEvaluatorOperator.




**Warning:**

The derivatives cannot be computed outside of the domain.




**Parameters:**


* `patches_deriv_1` A [**MultipatchType**](classMultipatchType.md) of DField to store the derivatives of the splines at the given coordinates. 
* `patches_coords` A [**MultipatchType**](classMultipatchType.md) of Field of Coordinate storing the coordinates of the meshes. 
* `patches_splines` A [**MultipatchType**](classMultipatchType.md) of DField storing the 2D spline coefficients. 




        

<hr>



### function deriv\_dim\_2 [1/2]

_Differentiate 2D splines (described by their spline coefficients) at a given coordinate along second dimension of interest._ 
```C++
template<class Coord>
inline KOKKOS_FUNCTION double MultipatchSplineEvaluator2D::deriv_dim_2 (
    Coord const & coord_eval,
    MultipatchSplineCoeff const & patches_splines
) const
```



See MultipatchSplineEvaluatorOperator.




**Warning:**

The derivative cannot be computed outside of the domain.




**Template parameters:**


* `StoringPatch` [**Patch**](structPatch.md) type where the given coordinate is stored. It does not mean that the coordinate is physically located on the patch. 



**Parameters:**


* `coord_eval` The coordinate where the spline is differentiated. 
* `patches_splines` A [**MultipatchType**](classMultipatchType.md) of DField storing the 2D spline coefficients.



**Returns:**

The derivative of the spline at the desired coordinate on the right patch. 





        

<hr>



### function deriv\_dim\_2 [2/2]

_Differentiate 2D splines (described by their spline coefficients) on a meshes along second dimension of interest._ 
```C++
inline void MultipatchSplineEvaluator2D::deriv_dim_2 (
    MultipatchValues const & patches_deriv_2,
    MultipatchCoordField const & patches_coords,
    MultipatchSplineCoeff const & patches_splines
) const
```



See MultipatchSplineEvaluatorOperator.




**Warning:**

The derivatives cannot be computed outside of the domain.




**Parameters:**


* `patches_deriv_2` A [**MultipatchType**](classMultipatchType.md) of DField to store the derivatives of the splines at the given coordinates. 
* `patches_coords` A [**MultipatchType**](classMultipatchType.md) of Field of Coordinate storing the coordinates of the meshes. 
* `patches_splines` A [**MultipatchType**](classMultipatchType.md) of DField storing the 2D spline coefficients. 




        

<hr>



### function integrate 

_Integration of splines (described by their spline coefficients)._ 
```C++
inline void MultipatchSplineEvaluator2D::integrate (
    DKokkosView_h< n_patches > const & integrals,
    MultipatchSplineCoeff const & patches_splines
) const
```



See MultipatchSplineEvaluatorOperator.




**Parameters:**


* `integrals` An Kokkos::View on host containing the integrals of each spline on each patch. 
* `patches_splines` A [**MultipatchType**](classMultipatchType.md) of DField storing the 2D spline coefficients. 




        

<hr>



### function operator() 

_Evaluate the 2D splines (described by their spline coefficients) at a given coordinate._ 
```C++
template<class Coord>
inline KOKKOS_FUNCTION double MultipatchSplineEvaluator2D::operator() (
    Coord const coord_eval,
    MultipatchSplineCoeff const & patches_splines
) const
```



The spline coefficients represent a 2D spline function defined on a B-splines (basis splines). They can be obtained via various methods, such as using a SplineBuilder2D or [**MultipatchSplineBuilder2D**](classMultipatchSplineBuilder2D.md). The coordinate is defined on one patch. But it can be physically located on another patch.


 

**Template parameters:**


* `StoringPatch` [**Patch**](structPatch.md) type where the given coordinate is stored. It does not mean that the coordinate is physically located on the patch. 



**Parameters:**


* `coord_eval` The coordinate where the spline is evaluated. 
* `patches_splines` A [**MultipatchType**](classMultipatchType.md) of DField storing the 2D spline coefficients. 



**Returns:**

The value of the spline at the desired coordinate on the right patch. 





        

<hr>



### function operator() 

_Evaluate 2D splines (described by their spline coefficients) on meshes._ 
```C++
inline void MultipatchSplineEvaluator2D::operator() (
    MultipatchValues const & patches_values,
    MultipatchCoordField const & patches_coords,
    MultipatchSplineCoeff const & patches_splines
) const
```



See MultipatchSplineEvaluatorOperator.




**Parameters:**


* `patches_values` A [**MultipatchType**](classMultipatchType.md) of DField to store the values of the splines at the given coordinates. 
* `patches_coords` A [**MultipatchType**](classMultipatchType.md) of Field of Coordinate storing the coordinates of the meshes. 
* `patches_splines` A [**MultipatchType**](classMultipatchType.md) of DField storing the 2D spline coefficients. 




        

<hr>



### function ~MultipatchSplineEvaluator2D 

```C++
MultipatchSplineEvaluator2D::~MultipatchSplineEvaluator2D () = default
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/spline/multipatch_spline_evaluator_2d.hpp`

