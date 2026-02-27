

# Struct InterpolationEvaluatorTraits

**template &lt;class Evaluator&gt;**



[**ClassList**](annotated.md) **>** [**InterpolationEvaluatorTraits**](structInterpolationEvaluatorTraits.md)



_A traits struct for accessing type aliases of an interpolation evaluator._ [More...](#detailed-description)

* `#include <i_interpolation_evaluator.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef typename Evaluator::template [**batched\_coeff\_idx\_range\_type**](structInterpolationEvaluatorTraits.md#typedef-batched_coeff_idx_range_type)&lt; D &gt; | [**batched\_coeff\_idx\_range\_type**](#typedef-batched_coeff_idx_range_type)  <br>_Batched domain with the evaluation grid replaced by coeff\_grid\_type._  |
| typedef typename Evaluator::template [**batched\_evaluation\_idx\_range\_type**](structInterpolationEvaluatorTraits.md#typedef-batched_evaluation_idx_range_type)&lt; BatchedInterpolationIdxRange &gt; | [**batched\_evaluation\_idx\_range\_type**](#typedef-batched_evaluation_idx_range_type)  <br>_Batched index range for the evaluation._  |
| typedef typename Evaluator::coeff\_grid\_type | [**coeff\_grid\_type**](#typedef-coeff_grid_type)  <br>_The discrete dimension for the interpolation coefficients._  |
| typedef typename Evaluator::data\_type | [**data\_type**](#typedef-data_type)  <br>_The data type that the data is saved on._  |
| typedef typename Evaluator::evaluation\_idx\_range\_type | [**evaluation\_idx\_range\_type**](#typedef-evaluation_idx_range_type)  <br>_The 1D index range for the evaluation mesh._  |
















































## Detailed Description


The primary template delegates to the evaluator's own type aliases. It handles the common case where an evaluator already uses the convention names (e.g. [**LagrangeEvaluator**](classLagrangeEvaluator.md) exposes `evaluation_idx_range_type` and `batched_coeff_idx_range_type` directly).


Specialise this struct to adapt external evaluators whose alias names differ (e.g. ddc::SplineEvaluator).




**Template parameters:**


* `Evaluator` The interpolation evaluator type. 




    
## Public Types Documentation




### typedef batched\_coeff\_idx\_range\_type 

_Batched domain with the evaluation grid replaced by coeff\_grid\_type._ 
```C++
using InterpolationEvaluatorTraits< Evaluator >::batched_coeff_idx_range_type =  typename Evaluator::template batched_coeff_idx_range_type<D>;
```




<hr>



### typedef batched\_evaluation\_idx\_range\_type 

_Batched index range for the evaluation._ 
```C++
using InterpolationEvaluatorTraits< Evaluator >::batched_evaluation_idx_range_type =  typename Evaluator::template batched_evaluation_idx_range_type< BatchedInterpolationIdxRange>;
```




<hr>



### typedef coeff\_grid\_type 

_The discrete dimension for the interpolation coefficients._ 
```C++
using InterpolationEvaluatorTraits< Evaluator >::coeff_grid_type =  typename Evaluator::coeff_grid_type;
```




<hr>



### typedef data\_type 

_The data type that the data is saved on._ 
```C++
using InterpolationEvaluatorTraits< Evaluator >::data_type =  typename Evaluator::data_type;
```




<hr>



### typedef evaluation\_idx\_range\_type 

_The 1D index range for the evaluation mesh._ 
```C++
using InterpolationEvaluatorTraits< Evaluator >::evaluation_idx_range_type =  typename Evaluator::evaluation_idx_range_type;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/interpolation/i_interpolation_evaluator.hpp`

