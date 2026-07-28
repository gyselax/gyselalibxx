

# Class LagrangeInterpolator

**template &lt;class ExecSpace, class Basis, class InterpGrid, class MinExtrapRule, class MaxExtrapRule, class DataType&gt;**



[**ClassList**](annotated.md) **>** [**LagrangeInterpolator**](classLagrangeInterpolator.md)



_An owning interpolation object that bundles a Lagrange builder and evaluator._ [More...](#detailed-description)

* `#include <lagrange_interpolation.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef [**IdentityInterpolationBuilder**](classIdentityInterpolationBuilder.md)&lt; ExecSpace, typename ExecSpace::memory\_space, DataType, InterpGrid, Basis &gt; | [**BuilderType**](#typedef-buildertype)  <br>_The_ [_**IdentityInterpolationBuilder**_](classIdentityInterpolationBuilder.md) _type built from the template parameters._ |
| typedef typename [**BuilderType::basis\_domain\_type**](classIdentityInterpolationBuilder.md#typedef-basis_domain_type) | [**CoeffGridType**](#typedef-coeffgridtype)  <br>_The discrete grid type used for the Lagrange coefficients (the Lagrange basis grid)._  |
| typedef [**LagrangeEvaluator**](classLagrangeEvaluator.md)&lt; ExecSpace, typename ExecSpace::memory\_space, DataType, Basis, InterpGrid, MinExtrapolationRule, MaxExtrapolationRule &gt; | [**EvaluatorType**](#typedef-evaluatortype)  <br>_The_ [_**LagrangeEvaluator**_](classLagrangeEvaluator.md) _type built from the template parameters._ |




















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**LagrangeInterpolator**](#function-lagrangeinterpolator-12) (IdxRange&lt; InterpGrid &gt; idx\_range=IdxRange&lt; InterpGrid &gt; {}) <br>_Construct a_ [_**LagrangeInterpolator**_](classLagrangeInterpolator.md) _._ |
|   | [**LagrangeInterpolator**](#function-lagrangeinterpolator-22) (IdxRange&lt; InterpGrid &gt; idx\_range, MinExtrapolationRule min\_extrapolation\_rule, MaxExtrapolationRule max\_extrapolation\_rule) <br>_Construct a_ [_**LagrangeInterpolator**_](classLagrangeInterpolator.md) _, specifying the extrapolation rules explicitly._ |
|  [**BuilderType**](classLagrangeInterpolator.md#typedef-buildertype) const & | [**get\_builder**](#function-get_builder) () const<br>_Return a const reference to the owned identity builder._  |
|  [**EvaluatorType**](classLagrangeInterpolator.md#typedef-evaluatortype) const & | [**get\_evaluator**](#function-get_evaluator) () const<br>_Return a const reference to the owned Lagrange evaluator._  |


## Public Static Functions

| Type | Name |
| ---: | :--- |
|  constexpr std::size\_t | [**rank**](#function-rank) () <br>_The number of interpolation dimensions._  |


























## Detailed Description


[**LagrangeInterpolator**](classLagrangeInterpolator.md) constructs and owns a matching [**IdentityInterpolationBuilder**](classIdentityInterpolationBuilder.md) and [**LagrangeEvaluator**](classLagrangeEvaluator.md) for a given dimension. It satisfies the concepts::Interpolation concept and is the recommended way to create a Lagrange interpolation for use with advection operators and similar algorithms.


The builder is an identity operation: it passes function values on the interpolation grid directly as coefficients to the evaluator, which then performs local polynomial reconstruction via the Lagrange basis.




**Template parameters:**


* `ExecSpace` The Kokkos execution space used for computations. 
* `Basis` The Lagrange basis type (uniform or non-uniform). 
* `InterpGrid` The discrete grid on which function values are provided. 
* `MinExtrapRule` The extrapolation rule applied below the lower boundary. This may be one of the tags in the [**ExtrapolationRule**](namespaceExtrapolationRule.md) namespace (e.g. [**ExtrapolationRule::Periodic**](structExtrapolationRule_1_1Periodic.md)) or a custom, already-concrete extrapolation rule class. 
* `MaxExtrapRule` The extrapolation rule applied above the upper boundary. See MinExtrapRule for the accepted forms. 
* `DataType` The floating-point type of the function values (default: double). 




    
## Public Types Documentation




### typedef BuilderType 

_The_ [_**IdentityInterpolationBuilder**_](classIdentityInterpolationBuilder.md) _type built from the template parameters._
```C++
using LagrangeInterpolator< ExecSpace, Basis, InterpGrid, MinExtrapRule, MaxExtrapRule, DataType >::BuilderType =  IdentityInterpolationBuilder< ExecSpace, typename ExecSpace::memory_space, DataType, InterpGrid, Basis>;
```




<hr>



### typedef CoeffGridType 

_The discrete grid type used for the Lagrange coefficients (the Lagrange basis grid)._ 
```C++
using LagrangeInterpolator< ExecSpace, Basis, InterpGrid, MinExtrapRule, MaxExtrapRule, DataType >::CoeffGridType =  typename BuilderType::basis_domain_type;
```




<hr>



### typedef EvaluatorType 

_The_ [_**LagrangeEvaluator**_](classLagrangeEvaluator.md) _type built from the template parameters._
```C++
using LagrangeInterpolator< ExecSpace, Basis, InterpGrid, MinExtrapRule, MaxExtrapRule, DataType >::EvaluatorType =  LagrangeEvaluator< ExecSpace, typename ExecSpace::memory_space, DataType, Basis, InterpGrid, MinExtrapolationRule, MaxExtrapolationRule>;
```




<hr>
## Public Functions Documentation




### function LagrangeInterpolator [1/2]

_Construct a_ [_**LagrangeInterpolator**_](classLagrangeInterpolator.md) _._
```C++
inline explicit LagrangeInterpolator::LagrangeInterpolator (
    IdxRange< InterpGrid > idx_range=IdxRange< InterpGrid > {}
) 
```



The extrapolation rules are initialised from the discrete space of `Basis`, so the corresponding ddc discrete space must be initialised before construction. No index range is required because the identity builder needs none. This overload is only available when both extrapolation rules can be built automatically (they are default-constructible, or the tag [**ExtrapolationRule::Constant**](structExtrapolationRule_1_1Constant.md) is used) - otherwise use the overload that takes the extrapolation rules explicitly.




**Parameters:**


* `idx_range` The index range on which the interpolator will act. This is unused but is included to match the [**SplineInterpolator**](classSplineInterpolator.md) interface. 




        

<hr>



### function LagrangeInterpolator [2/2]

_Construct a_ [_**LagrangeInterpolator**_](classLagrangeInterpolator.md) _, specifying the extrapolation rules explicitly._
```C++
inline explicit LagrangeInterpolator::LagrangeInterpolator (
    IdxRange< InterpGrid > idx_range,
    MinExtrapolationRule min_extrapolation_rule,
    MaxExtrapolationRule max_extrapolation_rule
) 
```



Use this overload when the chosen extrapolation rule cannot be built automatically, e.g. a custom extrapolation rule that is not default-constructible and is not [**ExtrapolationRule::Constant**](structExtrapolationRule_1_1Constant.md).




**Parameters:**


* `idx_range` The index range on which the interpolator will act. This is unused but is included to match the [**SplineInterpolator**](classSplineInterpolator.md) interface. 
* `min_extrapolation_rule` The extrapolation rule to use below the lower boundary. 
* `max_extrapolation_rule` The extrapolation rule to use above the upper boundary. 




        

<hr>



### function get\_builder 

_Return a const reference to the owned identity builder._ 
```C++
inline BuilderType const & LagrangeInterpolator::get_builder () const
```





**Returns:**

The BuilderType instance. 





        

<hr>



### function get\_evaluator 

_Return a const reference to the owned Lagrange evaluator._ 
```C++
inline EvaluatorType const & LagrangeInterpolator::get_evaluator () const
```





**Returns:**

The EvaluatorType instance. 





        

<hr>
## Public Static Functions Documentation




### function rank 

_The number of interpolation dimensions._ 
```C++
static inline constexpr std::size_t LagrangeInterpolator::rank () 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/interpolation/lagrange_interpolation.hpp`

