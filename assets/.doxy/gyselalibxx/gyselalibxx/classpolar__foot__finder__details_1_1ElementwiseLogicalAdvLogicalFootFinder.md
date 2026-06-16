

# Class polar\_foot\_finder\_details::ElementwiseLogicalAdvLogicalFootFinder

**template &lt;class [**GridR**](structGridR.md), class [**GridTheta**](structGridTheta.md), class IdxRangeOperator, class RThetaAdvectionEvaluator, class AdvecCoefField, class TimeStepper&gt;**



[**ClassList**](annotated.md) **>** [**polar\_foot\_finder\_details**](namespacepolar__foot__finder__details.md) **>** [**ElementwiseLogicalAdvLogicalFootFinder**](classpolar__foot__finder__details_1_1ElementwiseLogicalAdvLogicalFootFinder.md)



_GPU-callable functor that finds characteristic feet for logical-space advection with foot-finding in logical_ \((r, \theta)\) _space._[More...](#detailed-description)

* `#include <logical_advection_logical_foot_finder.hpp>`





































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**ElementwiseLogicalAdvLogicalFootFinder**](#function-elementwiselogicaladvlogicalfootfinder) (RThetaAdvectionEvaluator const & evaluator\_advection\_field, TimeStepper const & time\_stepper, AdvecCoefField const & advection\_field\_coefs, IdxRange&lt; [**GridTheta**](structGridTheta.md) &gt; idx\_range\_theta, double dt) <br>_Construct an_ [_**ElementwiseLogicalAdvLogicalFootFinder**_](classpolar__foot__finder__details_1_1ElementwiseLogicalAdvLogicalFootFinder.md) _._ |
|  KOKKOS\_FUNCTION CoordRTheta | [**operator()**](#function-operator) (IdxOperator const idx) const<br>_Find the foot of the characteristic at a single grid index._  |




























## Detailed Description


Constructed by [**ElementwiseLogicalAdvLogicalFootFinderMem**](classpolar__foot__finder__details_1_1ElementwiseLogicalAdvLogicalFootFinderMem.md) and intended to be copied to the device and called elementwise inside a Kokkos kernel. 


    
## Public Functions Documentation




### function ElementwiseLogicalAdvLogicalFootFinder 

_Construct an_ [_**ElementwiseLogicalAdvLogicalFootFinder**_](classpolar__foot__finder__details_1_1ElementwiseLogicalAdvLogicalFootFinder.md) _._
```C++
inline polar_foot_finder_details::ElementwiseLogicalAdvLogicalFootFinder::ElementwiseLogicalAdvLogicalFootFinder (
    RThetaAdvectionEvaluator const & evaluator_advection_field,
    TimeStepper const & time_stepper,
    AdvecCoefField const & advection_field_coefs,
    IdxRange< GridTheta > idx_range_theta,
    double dt
) 
```





**Parameters:**


* `evaluator_advection_field` The evaluator for the spline representation of the advection field. 
* `time_stepper` The time integration method. 
* `advection_field_coefs` A view of the spline coefficients of the advection field. 
* `idx_range_theta` The poloidal index range, used to wrap the angular coordinate into the periodic domain after each time step. 
* `dt` The time step for the characteristic equation. 




        

<hr>



### function operator() 

_Find the foot of the characteristic at a single grid index._ 
```C++
inline KOKKOS_FUNCTION CoordRTheta polar_foot_finder_details::ElementwiseLogicalAdvLogicalFootFinder::operator() (
    IdxOperator const idx
) const
```



Solves the characteristic equation over \(dt\) using the stored time stepper and returns the resulting \((r, \theta)\) coordinate.




**Parameters:**


* `idx` The operator index, encoding both batch and \((r, \theta)\) indices.



**Returns:**

The \((r, \theta)\) coordinate of the characteristic foot. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/advection/polar_foot_finders/logical_advection_logical_foot_finder.hpp`

