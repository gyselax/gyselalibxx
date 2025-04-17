

# Class PoissonLikeRHSFunction

**template &lt;class RadialExtrapolationRule&gt;**



[**ClassList**](annotated.md) **>** [**PoissonLikeRHSFunction**](classPoissonLikeRHSFunction.md)



_Type of right-hand side (rhs) function of the Poisson equation._ [More...](#detailed-description)

* `#include <poisson_like_rhs_function.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef ddc::SplineEvaluator2D&lt; Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace, [**BSplinesR**](structBSplinesR.md), [**BSplinesTheta**](structBSplinesTheta.md), [**GridR**](structGridR.md), [**GridTheta**](structGridTheta.md), RadialExtrapolationRule, RadialExtrapolationRule, ddc::PeriodicExtrapolationRule&lt; [**Theta**](structTheta.md) &gt;, ddc::PeriodicExtrapolationRule&lt; [**Theta**](structTheta.md) &gt; &gt; | [**evaluator\_type**](#typedef-evaluator_type)  <br>_The type of the 2D Spline Evaluator used by this class._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**PoissonLikeRHSFunction**](#function-poissonlikerhsfunction) (host\_t&lt; ConstSpline2D &gt; coefs, [**evaluator\_type**](classPoissonLikeRHSFunction.md#typedef-evaluator_type) const & evaluator) <br>_Instantiate a_ [_**PoissonLikeRHSFunction**_](classPoissonLikeRHSFunction.md) _._ |
|  double | [**operator()**](#function-operator) (CoordRTheta const & coord\_rtheta) const<br>_Get the value of the function at a given coordinate._  |




























## Detailed Description




**Template parameters:**


* `RadialExtrapolationRule` The extrapolation rule applied at the outer radial bound. 




    
## Public Types Documentation




### typedef evaluator\_type 

_The type of the 2D Spline Evaluator used by this class._ 
```C++
using PoissonLikeRHSFunction< RadialExtrapolationRule >::evaluator_type =  ddc::SplineEvaluator2D< Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace, BSplinesR, BSplinesTheta, GridR, GridTheta, RadialExtrapolationRule, RadialExtrapolationRule, ddc::PeriodicExtrapolationRule<Theta>, ddc::PeriodicExtrapolationRule<Theta> >;
```




<hr>
## Public Functions Documentation




### function PoissonLikeRHSFunction 

_Instantiate a_ [_**PoissonLikeRHSFunction**_](classPoissonLikeRHSFunction.md) _._
```C++
inline PoissonLikeRHSFunction::PoissonLikeRHSFunction (
    host_t< ConstSpline2D > coefs,
    evaluator_type const & evaluator
) 
```





**Parameters:**


* `coefs` The B-splines coefficients of the right-hand side function. 
* `evaluator` Evaluator on B-splines. 




        

<hr>



### function operator() 

_Get the value of the function at a given coordinate._ 
```C++
inline double PoissonLikeRHSFunction::operator() (
    CoordRTheta const & coord_rtheta
) const
```





**Parameters:**


* `coord_rtheta` Polar coordinate where we want to evaluate the rhs function.



**Returns:**

A double with the value of the rhs at the given coordinate. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/geometryRTheta/poisson/poisson_like_rhs_function.hpp`

