

# Class SplineInterpolator

**template &lt;class ExecSpace, class Basis, class InterpGrid, class MinExtrapRule, class MaxExtrapRule, ddc::SplineBuilderClosure MinBound, ddc::SplineBuilderClosure MaxBound, ddc::SplineSolver Solver&gt;**



[**ClassList**](annotated.md) **>** [**SplineInterpolator**](classSplineInterpolator.md)



_An owning interpolation object that bundles a spline builder and evaluator._ [More...](#detailed-description)

* `#include <spline_interpolation.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef ddc::SplineBuilder&lt; ExecSpace, typename ExecSpace::memory\_space, Basis, InterpGrid, MinBound, MaxBound, Solver &gt; | [**BuilderType**](#typedef-buildertype)  <br>_The ddc::SplineBuilder type built from the template parameters._  |
| typedef ddc::SplineEvaluator&lt; ExecSpace, typename ExecSpace::memory\_space, Basis, InterpGrid, MinExtrapolationRule, MaxExtrapolationRule &gt; | [**EvaluatorType**](#typedef-evaluatortype)  <br>_The ddc::SplineEvaluator type built from the template parameters._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**SplineInterpolator**](#function-splineinterpolator-14) (std::string const & label, IdxRange&lt; InterpGrid &gt; idx\_range) <br>_Construct a_ [_**SplineInterpolator**_](classSplineInterpolator.md) _on the given interpolation index range._ |
|   | [**SplineInterpolator**](#function-splineinterpolator-24) (IdxRange&lt; InterpGrid &gt; idx\_range) <br>_Construct a_ [_**SplineInterpolator**_](classSplineInterpolator.md) _on the given interpolation index range._ |
|   | [**SplineInterpolator**](#function-splineinterpolator-34) (std::string const & label, IdxRange&lt; InterpGrid &gt; idx\_range, MinExtrapolationRule min\_extrapolation\_rule, MaxExtrapolationRule max\_extrapolation\_rule) <br>_Construct a_ [_**SplineInterpolator**_](classSplineInterpolator.md) _on the given interpolation index range, specifying the extrapolation rules explicitly._ |
|   | [**SplineInterpolator**](#function-splineinterpolator-44) (IdxRange&lt; InterpGrid &gt; idx\_range, MinExtrapolationRule min\_extrapolation\_rule, MaxExtrapolationRule max\_extrapolation\_rule) <br>_Construct a_ [_**SplineInterpolator**_](classSplineInterpolator.md) _on the given interpolation index range, specifying the extrapolation rules explicitly._ |
|  [**BuilderType**](classSplineInterpolator.md#typedef-buildertype) const & | [**get\_builder**](#function-get_builder) () const<br>_Return a const reference to the owned spline builder._  |
|  [**EvaluatorType**](classSplineInterpolator.md#typedef-evaluatortype) const & | [**get\_evaluator**](#function-get_evaluator) () const<br>_Return a const reference to the owned spline evaluator._  |


## Public Static Functions

| Type | Name |
| ---: | :--- |
|  constexpr std::size\_t | [**rank**](#function-rank) () <br>_The number of interpolation dimensions._  |


























## Detailed Description


[**SplineInterpolator**](classSplineInterpolator.md) constructs and owns a matching ddc::SplineBuilder and ddc::SplineEvaluator for a given dimension. It satisfies the concepts::Interpolation concept and is the recommended way to create a spline interpolation for use with advection operators and similar algorithms.


The boundary condition (MinBound / MaxBound) and extrapolation rule (MinExtrapRule / MaxExtrapRule) must be consistent: both must be PERIODIC for periodic dimensions and both must be non-PERIODIC for non-periodic dimensions.




**Template parameters:**


* `ExecSpace` The Kokkos execution space used for computations. 
* `Basis` The B-spline basis type (uniform or non-uniform). 
* `InterpGrid` The discrete grid on which function values are provided. 
* `MinExtrapRule` The extrapolation rule applied below the lower boundary. This may be one of the tags in the [**ExtrapolationRule**](namespaceExtrapolationRule.md) namespace (e.g. [**ExtrapolationRule::Periodic**](structExtrapolationRule_1_1Periodic.md)) or a custom, already-concrete extrapolation rule class. 
* `MaxExtrapRule` The extrapolation rule applied above the upper boundary. See MinExtrapRule for the accepted forms. 
* `MinBound` The ddc::SplineBuilderClosure at the lower boundary of the spline builder. 
* `MaxBound` The ddc::SplineBuilderClosure at the upper boundary of the spline builder. 
* `Solver` The spline solver backend (default: LAPACK). 




    
## Public Types Documentation




### typedef BuilderType 

_The ddc::SplineBuilder type built from the template parameters._ 
```C++
using SplineInterpolator< ExecSpace, Basis, InterpGrid, MinExtrapRule, MaxExtrapRule, MinBound, MaxBound, Solver >::BuilderType =  ddc::SplineBuilder< ExecSpace, typename ExecSpace::memory_space, Basis, InterpGrid, MinBound, MaxBound, Solver>;
```




<hr>



### typedef EvaluatorType 

_The ddc::SplineEvaluator type built from the template parameters._ 
```C++
using SplineInterpolator< ExecSpace, Basis, InterpGrid, MinExtrapRule, MaxExtrapRule, MinBound, MaxBound, Solver >::EvaluatorType =  ddc::SplineEvaluator< ExecSpace, typename ExecSpace::memory_space, Basis, InterpGrid, MinExtrapolationRule, MaxExtrapolationRule>;
```




<hr>
## Public Functions Documentation




### function SplineInterpolator [1/4]

_Construct a_ [_**SplineInterpolator**_](classSplineInterpolator.md) _on the given interpolation index range._
```C++
inline explicit SplineInterpolator::SplineInterpolator (
    std::string const & label,
    IdxRange< InterpGrid > idx_range
) 
```



The extrapolation rules are initialised from the discrete space of `Basis`, so the corresponding ddc discrete space must be initialised before construction. This overload is only available when both extrapolation rules can be built automatically (they are default-constructible, or the tag [**ExtrapolationRule::Constant**](structExtrapolationRule_1_1Constant.md) is used) - otherwise use the overload that takes the extrapolation rules explicitly.




**Parameters:**


* `label` A label used to tag parallel regions and memory allocations for profiling. 
* `idx_range` The 1D interpolation index range passed to the builder. 




        

<hr>



### function SplineInterpolator [2/4]

_Construct a_ [_**SplineInterpolator**_](classSplineInterpolator.md) _on the given interpolation index range._
```C++
inline explicit SplineInterpolator::SplineInterpolator (
    IdxRange< InterpGrid > idx_range
) 
```



The extrapolation rules are initialised from the discrete space of `Basis`, so the corresponding ddc discrete space must be initialised before construction. This overload is only available when both extrapolation rules can be built automatically (they are default-constructible, or the tag [**ExtrapolationRule::Constant**](structExtrapolationRule_1_1Constant.md) is used) - otherwise use the overload that takes the extrapolation rules explicitly.




**Parameters:**


* `idx_range` The 1D interpolation index range passed to the builder. 




        

<hr>



### function SplineInterpolator [3/4]

_Construct a_ [_**SplineInterpolator**_](classSplineInterpolator.md) _on the given interpolation index range, specifying the extrapolation rules explicitly._
```C++
inline explicit SplineInterpolator::SplineInterpolator (
    std::string const & label,
    IdxRange< InterpGrid > idx_range,
    MinExtrapolationRule min_extrapolation_rule,
    MaxExtrapolationRule max_extrapolation_rule
) 
```



Use this overload when the chosen extrapolation rule cannot be built automatically, e.g. a custom extrapolation rule that is not default-constructible and is not [**ExtrapolationRule::Constant**](structExtrapolationRule_1_1Constant.md).




**Parameters:**


* `label` A label used to tag parallel regions and memory allocations for profiling. 
* `idx_range` The 1D interpolation index range passed to the builder. 
* `min_extrapolation_rule` The extrapolation rule to use below the lower boundary. 
* `max_extrapolation_rule` The extrapolation rule to use above the upper boundary. 




        

<hr>



### function SplineInterpolator [4/4]

_Construct a_ [_**SplineInterpolator**_](classSplineInterpolator.md) _on the given interpolation index range, specifying the extrapolation rules explicitly._
```C++
inline explicit SplineInterpolator::SplineInterpolator (
    IdxRange< InterpGrid > idx_range,
    MinExtrapolationRule min_extrapolation_rule,
    MaxExtrapolationRule max_extrapolation_rule
) 
```



Use this overload when the chosen extrapolation rule cannot be built automatically, e.g. a custom extrapolation rule that is not default-constructible and is not [**ExtrapolationRule::Constant**](structExtrapolationRule_1_1Constant.md).




**Parameters:**


* `idx_range` The 1D interpolation index range passed to the builder. 
* `min_extrapolation_rule` The extrapolation rule to use below the lower boundary. 
* `max_extrapolation_rule` The extrapolation rule to use above the upper boundary. 




        

<hr>



### function get\_builder 

_Return a const reference to the owned spline builder._ 
```C++
inline BuilderType const & SplineInterpolator::get_builder () const
```





**Returns:**

The BuilderType instance. 





        

<hr>



### function get\_evaluator 

_Return a const reference to the owned spline evaluator._ 
```C++
inline EvaluatorType const & SplineInterpolator::get_evaluator () const
```





**Returns:**

The EvaluatorType instance. 





        

<hr>
## Public Static Functions Documentation




### function rank 

_The number of interpolation dimensions._ 
```C++
static inline constexpr std::size_t SplineInterpolator::rank () 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/interpolation/spline_interpolation.hpp`

