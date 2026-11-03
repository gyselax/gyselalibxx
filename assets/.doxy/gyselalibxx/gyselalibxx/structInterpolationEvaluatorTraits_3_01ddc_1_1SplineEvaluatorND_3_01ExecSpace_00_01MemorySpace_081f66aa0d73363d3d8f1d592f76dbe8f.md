

# Struct InterpolationEvaluatorTraits&lt; ddc::SplineEvaluatorND&lt; ExecSpace, MemorySpace, ddc::TypeSeq&lt; BSplines... &gt;, ddc::TypeSeq&lt; EvaluationDDim... &gt;, ddc::TypeSeq&lt; ExtrapolationRule... &gt; &gt; &gt;

**template &lt;class ExecSpace, class MemorySpace, class... BSplines, class... EvaluationDDim, class... ExtrapolationRule&gt;**



[**ClassList**](annotated.md) **>** [**InterpolationEvaluatorTraits&lt; ddc::SplineEvaluatorND&lt; ExecSpace, MemorySpace, ddc::TypeSeq&lt; BSplines... &gt;, ddc::TypeSeq&lt; EvaluationDDim... &gt;, ddc::TypeSeq&lt; ExtrapolationRule... &gt; &gt; &gt;**](structInterpolationEvaluatorTraits_3_01ddc_1_1SplineEvaluatorND_3_01ExecSpace_00_01MemorySpace_081f66aa0d73363d3d8f1d592f76dbe8f.md)



_Specialisation of_ [_**InterpolationEvaluatorTraits**_](structInterpolationEvaluatorTraits.md) _for ddc::SplineEvaluatorND._[More...](#detailed-description)

* `#include <i_interpolation_evaluator.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef typename Evaluator::template batched\_spline\_domain\_type&lt; BatchedInterpolationIdxRange &gt; | [**batched\_coeff\_idx\_range\_type**](#typedef-batched_coeff_idx_range_type)  <br>_Batched domain with the evaluation grid replaced by BSplines._  |
| typedef typename Evaluator::template batched\_evaluation\_domain\_type&lt; BatchedInterpolationIdxRange &gt; | [**batched\_evaluation\_idx\_range\_type**](#typedef-batched_evaluation_idx_range_type)  <br>_Batched index range for the evaluation._  |
| typedef decltype(make\_coeff\_idx\_range\_type(std::make\_index\_sequence&lt; sizeof...(BSplines)&gt;())) | [**coeff\_idx\_range\_type**](#typedef-coeff_idx_range_type)  <br>_The type of the ND index range on which the interpolation coefficients are defined._  |
| typedef Coord&lt; typename EvaluationDDim::continuous\_dimension\_type... &gt; | [**coord\_type**](#typedef-coord_type)  <br>_The ND coordinate type corresponding to the evaluation mesh._  |
| typedef double | [**data\_type**](#typedef-data_type)  <br>_The data type that the data is saved on._  |
| typedef decltype(make\_evaluation\_idx\_range\_type( std::make\_index\_sequence&lt; sizeof...(BSplines)&gt;())) | [**evaluation\_idx\_range\_type**](#typedef-evaluation_idx_range_type)  <br>_The ND index range for the evaluation mesh._  |






















## Public Static Functions

| Type | Name |
| ---: | :--- |
|  std::size\_t | [**rank**](#function-rank) () <br>_The number of interpolation dimensions._  |


























## Detailed Description


ddc::SplineEvaluatorND uses different alias names from the InterpolationEvaluator convention. This specialisation provides the mapping so that ddc::SplineEvaluatorND can be used directly as an InterpolationEvaluator.


Mapping: evaluation\_domain\_type&lt;0, ...&gt; -&gt; [**evaluation\_idx\_range\_type**](structInterpolationEvaluatorTraits_3_01ddc_1_1SplineEvaluatorND_3_01ExecSpace_00_01MemorySpace_081f66aa0d73363d3d8f1d592f76dbe8f.md#typedef-evaluation_idx_range_type) spline\_domain\_type&lt;0, ...&gt; -&gt; [**coeff\_idx\_range\_type**](structInterpolationEvaluatorTraits_3_01ddc_1_1SplineEvaluatorND_3_01ExecSpace_00_01MemorySpace_081f66aa0d73363d3d8f1d592f76dbe8f.md#typedef-coeff_idx_range_type) batched\_spline\_domain\_type&lt;D&gt; -&gt; [**batched\_coeff\_idx\_range\_type&lt;D&gt;**](structInterpolationEvaluatorTraits_3_01ddc_1_1SplineEvaluatorND_3_01ExecSpace_00_01MemorySpace_081f66aa0d73363d3d8f1d592f76dbe8f.md#typedef-batched_coeff_idx_range_type) 


    
## Public Types Documentation




### typedef batched\_coeff\_idx\_range\_type 

_Batched domain with the evaluation grid replaced by BSplines._ 
```C++
using InterpolationEvaluatorTraits< ddc::SplineEvaluatorND< ExecSpace, MemorySpace, ddc::TypeSeq< BSplines... >, ddc::TypeSeq< EvaluationDDim... >, ddc::TypeSeq< ExtrapolationRule... > > >::batched_coeff_idx_range_type = 
            typename Evaluator::template batched_spline_domain_type<BatchedInterpolationIdxRange>;
```




<hr>



### typedef batched\_evaluation\_idx\_range\_type 

_Batched index range for the evaluation._ 
```C++
using InterpolationEvaluatorTraits< ddc::SplineEvaluatorND< ExecSpace, MemorySpace, ddc::TypeSeq< BSplines... >, ddc::TypeSeq< EvaluationDDim... >, ddc::TypeSeq< ExtrapolationRule... > > >::batched_evaluation_idx_range_type = 
            typename Evaluator::template batched_evaluation_domain_type<
                    BatchedInterpolationIdxRange>;
```




<hr>



### typedef coeff\_idx\_range\_type 

_The type of the ND index range on which the interpolation coefficients are defined._ 
```C++
using InterpolationEvaluatorTraits< ddc::SplineEvaluatorND< ExecSpace, MemorySpace, ddc::TypeSeq< BSplines... >, ddc::TypeSeq< EvaluationDDim... >, ddc::TypeSeq< ExtrapolationRule... > > >::coeff_idx_range_type =  decltype(make_coeff_idx_range_type(std::make_index_sequence<sizeof...(BSplines)>()));
```




<hr>



### typedef coord\_type 

_The ND coordinate type corresponding to the evaluation mesh._ 
```C++
using InterpolationEvaluatorTraits< ddc::SplineEvaluatorND< ExecSpace, MemorySpace, ddc::TypeSeq< BSplines... >, ddc::TypeSeq< EvaluationDDim... >, ddc::TypeSeq< ExtrapolationRule... > > >::coord_type =  Coord<typename EvaluationDDim::continuous_dimension_type...>;
```




<hr>



### typedef data\_type 

_The data type that the data is saved on._ 
```C++
using InterpolationEvaluatorTraits< ddc::SplineEvaluatorND< ExecSpace, MemorySpace, ddc::TypeSeq< BSplines... >, ddc::TypeSeq< EvaluationDDim... >, ddc::TypeSeq< ExtrapolationRule... > > >::data_type =  double;
```




<hr>



### typedef evaluation\_idx\_range\_type 

_The ND index range for the evaluation mesh._ 
```C++
using InterpolationEvaluatorTraits< ddc::SplineEvaluatorND< ExecSpace, MemorySpace, ddc::TypeSeq< BSplines... >, ddc::TypeSeq< EvaluationDDim... >, ddc::TypeSeq< ExtrapolationRule... > > >::evaluation_idx_range_type =  decltype(make_evaluation_idx_range_type(
            std::make_index_sequence<sizeof...(BSplines)>()));
```




<hr>
## Public Static Functions Documentation




### function rank 

_The number of interpolation dimensions._ 
```C++
static inline std::size_t InterpolationEvaluatorTraits< ddc::SplineEvaluatorND< ExecSpace, MemorySpace, ddc::TypeSeq< BSplines... >, ddc::TypeSeq< EvaluationDDim... >, ddc::TypeSeq< ExtrapolationRule... > > >::rank () 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/interpolation/i_interpolation_evaluator.hpp`

