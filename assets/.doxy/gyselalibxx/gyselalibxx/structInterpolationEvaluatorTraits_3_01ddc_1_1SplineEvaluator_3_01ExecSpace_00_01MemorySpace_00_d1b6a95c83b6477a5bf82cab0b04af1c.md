

# Struct InterpolationEvaluatorTraits&lt; ddc::SplineEvaluator&lt; ExecSpace, MemorySpace, BSplines, EvaluationDDim, LowerExtrapolationRule, UpperExtrapolationRule &gt; &gt;

**template &lt;class ExecSpace, class MemorySpace, class BSplines, class EvaluationDDim, class LowerExtrapolationRule, class UpperExtrapolationRule&gt;**



[**ClassList**](annotated.md) **>** [**InterpolationEvaluatorTraits&lt; ddc::SplineEvaluator&lt; ExecSpace, MemorySpace, BSplines, EvaluationDDim, LowerExtrapolationRule, UpperExtrapolationRule &gt; &gt;**](structInterpolationEvaluatorTraits_3_01ddc_1_1SplineEvaluator_3_01ExecSpace_00_01MemorySpace_00_d1b6a95c83b6477a5bf82cab0b04af1c.md)



_Specialisation of_ [_**InterpolationEvaluatorTraits**_](structInterpolationEvaluatorTraits.md) _for ddc::SplineEvaluator._[More...](#detailed-description)

* `#include <i_interpolation_evaluator.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef typename Evaluator::template batched\_spline\_domain\_type&lt; BatchedInterpolationIdxRange &gt; | [**batched\_coeff\_idx\_range\_type**](#typedef-batched_coeff_idx_range_type)  <br>_Batched domain with the evaluation grid replaced by BSplines._  |
| typedef typename Evaluator::template batched\_evaluation\_domain\_type&lt; BatchedInterpolationIdxRange &gt; | [**batched\_evaluation\_idx\_range\_type**](#typedef-batched_evaluation_idx_range_type)  <br>_Batched index range for the evaluation._  |
| typedef typename Evaluator::spline\_domain\_type | [**coeff\_idx\_range\_type**](#typedef-coeff_idx_range_type)  <br>_The type of the ND index range on which the interpolation coefficients are defined._  |
| typedef Coord&lt; typename EvaluationDDim::continuous\_dimension\_type &gt; | [**coord\_type**](#typedef-coord_type)  <br>_The 1D coordinate type corresponding to the evaluation mesh._  |
| typedef double | [**data\_type**](#typedef-data_type)  <br>_The data type that the data is saved on._  |
| typedef typename Evaluator::evaluation\_domain\_type | [**evaluation\_idx\_range\_type**](#typedef-evaluation_idx_range_type)  <br>_The 1D index range for the evaluation mesh._  |






















## Public Static Functions

| Type | Name |
| ---: | :--- |
|  constexpr std::size\_t | [**rank**](#function-rank) () <br>_The number of interpolation dimensions (always 1 for SplineEvaluator)._  |


























## Detailed Description


ddc::SplineEvaluator uses different alias names from the InterpolationEvaluator convention. This specialisation provides the mapping so that ddc::SplineEvaluator can be used directly as an InterpolationEvaluator.


Mapping: evaluation\_discrete\_dimension\_type -&gt; (defines evaluation\_idx\_range\_type) evaluation\_domain\_type -&gt; evaluation\_idx\_range\_type spline\_domain\_type -&gt; coeff\_idx\_range\_type batched\_spline\_domain\_type&lt;D&gt; -&gt; batched\_coeff\_idx\_range\_type&lt;D&gt; 


    
## Public Types Documentation




### typedef batched\_coeff\_idx\_range\_type 

_Batched domain with the evaluation grid replaced by BSplines._ 
```C++
using InterpolationEvaluatorTraits< ddc::SplineEvaluator< ExecSpace, MemorySpace, BSplines, EvaluationDDim, LowerExtrapolationRule, UpperExtrapolationRule > >::batched_coeff_idx_range_type =  typename Evaluator::template batched_spline_domain_type<BatchedInterpolationIdxRange>;
```




<hr>



### typedef batched\_evaluation\_idx\_range\_type 

_Batched index range for the evaluation._ 
```C++
using InterpolationEvaluatorTraits< ddc::SplineEvaluator< ExecSpace, MemorySpace, BSplines, EvaluationDDim, LowerExtrapolationRule, UpperExtrapolationRule > >::batched_evaluation_idx_range_type =  typename Evaluator::template batched_evaluation_domain_type< BatchedInterpolationIdxRange>;
```




<hr>



### typedef coeff\_idx\_range\_type 

_The type of the ND index range on which the interpolation coefficients are defined._ 
```C++
using InterpolationEvaluatorTraits< ddc::SplineEvaluator< ExecSpace, MemorySpace, BSplines, EvaluationDDim, LowerExtrapolationRule, UpperExtrapolationRule > >::coeff_idx_range_type =  typename Evaluator::spline_domain_type;
```




<hr>



### typedef coord\_type 

_The 1D coordinate type corresponding to the evaluation mesh._ 
```C++
using InterpolationEvaluatorTraits< ddc::SplineEvaluator< ExecSpace, MemorySpace, BSplines, EvaluationDDim, LowerExtrapolationRule, UpperExtrapolationRule > >::coord_type =  Coord<typename EvaluationDDim::continuous_dimension_type>;
```




<hr>



### typedef data\_type 

_The data type that the data is saved on._ 
```C++
using InterpolationEvaluatorTraits< ddc::SplineEvaluator< ExecSpace, MemorySpace, BSplines, EvaluationDDim, LowerExtrapolationRule, UpperExtrapolationRule > >::data_type =  double;
```




<hr>



### typedef evaluation\_idx\_range\_type 

_The 1D index range for the evaluation mesh._ 
```C++
using InterpolationEvaluatorTraits< ddc::SplineEvaluator< ExecSpace, MemorySpace, BSplines, EvaluationDDim, LowerExtrapolationRule, UpperExtrapolationRule > >::evaluation_idx_range_type =  typename Evaluator::evaluation_domain_type;
```




<hr>
## Public Static Functions Documentation




### function rank 

_The number of interpolation dimensions (always 1 for SplineEvaluator)._ 
```C++
static inline constexpr std::size_t InterpolationEvaluatorTraits< ddc::SplineEvaluator< ExecSpace, MemorySpace, BSplines, EvaluationDDim, LowerExtrapolationRule, UpperExtrapolationRule > >::rank () 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/interpolation/i_interpolation_evaluator.hpp`

