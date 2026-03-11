

# Class LagrangeInterpolator

**template &lt;class ExecSpace, class Basis, class InterpGrid, ExtrapolationRule MinExtrapRule, ExtrapolationRule MaxExtrapRule, ddc::BoundCond MinBound, ddc::BoundCond MaxBound, class DataType&gt;**



[**ClassList**](annotated.md) **>** [**LagrangeInterpolator**](classLagrangeInterpolator.md)



_An owning interpolation object that bundles a Lagrange builder and evaluator._ [More...](#detailed-description)

* `#include <lagrange_interpolation.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef [**IdentityInterpolationBuilder**](classIdentityInterpolationBuilder.md)&lt; ExecSpace, typename ExecSpace::memory\_space, DataType, InterpGrid, Basis &gt; | [**BuilderType**](#typedef-buildertype)  <br>_The_ [_**IdentityInterpolationBuilder**_](classIdentityInterpolationBuilder.md) _type built from the template parameters._ |
| typedef typename [**BuilderType::basis\_domain\_type**](classIdentityInterpolationBuilder.md#typedef-basis_domain_type) | [**CoeffGridType**](#typedef-coeffgridtype)  <br>_The discrete grid type used for the Lagrange coefficients (the Lagrange basis grid)._  |
| typedef [**LagrangeEvaluator**](classLagrangeEvaluator.md)&lt; ExecSpace, typename ExecSpace::memory\_space, DataType, Basis, InterpGrid, extrapolation\_rule\_t&lt; MinExtrapRule, [**CoeffGridType**](classLagrangeInterpolator.md#typedef-coeffgridtype) &gt;, extrapolation\_rule\_t&lt; MaxExtrapRule, [**CoeffGridType**](classLagrangeInterpolator.md#typedef-coeffgridtype) &gt; &gt; | [**EvaluatorType**](#typedef-evaluatortype)  <br>_The_ [_**LagrangeEvaluator**_](classLagrangeEvaluator.md) _type built from the template parameters._ |




















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**LagrangeInterpolator**](#function-lagrangeinterpolator) () <br>_Construct a_ [_**LagrangeInterpolator**_](classLagrangeInterpolator.md) _._ |
|  [**BuilderType**](classLagrangeInterpolator.md#typedef-buildertype) const & | [**get\_builder**](#function-get_builder) () const<br>_Return a const reference to the owned identity builder._  |
|  [**EvaluatorType**](classLagrangeInterpolator.md#typedef-evaluatortype) const & | [**get\_evaluator**](#function-get_evaluator) () const<br>_Return a const reference to the owned Lagrange evaluator._  |




























## Detailed Description


[**LagrangeInterpolator**](classLagrangeInterpolator.md) constructs and owns a matching [**IdentityInterpolationBuilder**](classIdentityInterpolationBuilder.md) and [**LagrangeEvaluator**](classLagrangeEvaluator.md) for a given dimension. It satisfies the concepts::Interpolation concept and is the recommended way to create a Lagrange interpolation for use with advection operators and similar algorithms.


The builder is an identity operation: it passes function values on the interpolation grid directly as coefficients to the evaluator, which then performs local polynomial reconstruction via the Lagrange basis.


The boundary condition (MinBound / MaxBound) and extrapolation rule (MinExtrapRule / MaxExtrapRule) must be consistent: both must be PERIODIC for periodic dimensions and both must be non-PERIODIC for non-periodic dimensions. Note: `CONSTANT` extrapolation is not supported for Lagrange interpolation.




**Template parameters:**


* `ExecSpace` The Kokkos execution space used for computations. 
* `Basis` The Lagrange basis type (uniform or non-uniform). 
* `InterpGrid` The discrete grid on which function values are provided. 
* `MinExtrapRule` The ExtrapolationRule applied below the lower boundary. 
* `MaxExtrapRule` The ExtrapolationRule applied above the upper boundary. 
* `MinBound` The ddc::BoundCond at the lower boundary (default: GREVILLE). This is included to have an interface interchangeable with SplineBuilder but is unused. 
* `MaxBound` The ddc::BoundCond at the upper boundary (default: GREVILLE). This is included to have an interface interchangeable with SplineBuilder but is unused. 
* `DataType` The floating-point type of the function values (default: double). 




    
## Public Types Documentation




### typedef BuilderType 

_The_ [_**IdentityInterpolationBuilder**_](classIdentityInterpolationBuilder.md) _type built from the template parameters._
```C++
using LagrangeInterpolator< ExecSpace, Basis, InterpGrid, MinExtrapRule, MaxExtrapRule, MinBound, MaxBound, DataType >::BuilderType =  IdentityInterpolationBuilder< ExecSpace, typename ExecSpace::memory_space, DataType, InterpGrid, Basis>;
```




<hr>



### typedef CoeffGridType 

_The discrete grid type used for the Lagrange coefficients (the Lagrange basis grid)._ 
```C++
using LagrangeInterpolator< ExecSpace, Basis, InterpGrid, MinExtrapRule, MaxExtrapRule, MinBound, MaxBound, DataType >::CoeffGridType =  typename BuilderType::basis_domain_type;
```




<hr>



### typedef EvaluatorType 

_The_ [_**LagrangeEvaluator**_](classLagrangeEvaluator.md) _type built from the template parameters._
```C++
using LagrangeInterpolator< ExecSpace, Basis, InterpGrid, MinExtrapRule, MaxExtrapRule, MinBound, MaxBound, DataType >::EvaluatorType =  LagrangeEvaluator< ExecSpace, typename ExecSpace::memory_space, DataType, Basis, InterpGrid, extrapolation_rule_t<MinExtrapRule, CoeffGridType>, extrapolation_rule_t<MaxExtrapRule, CoeffGridType> >;
```




<hr>
## Public Functions Documentation




### function LagrangeInterpolator 

_Construct a_ [_**LagrangeInterpolator**_](classLagrangeInterpolator.md) _._
```C++
inline LagrangeInterpolator::LagrangeInterpolator () 
```



The extrapolation rules are initialised from the discrete space of `Basis`, so the corresponding ddc discrete space must be initialised before construction. No index range is required because the identity builder needs none. 


        

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

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/interpolation/lagrange_interpolation.hpp`

