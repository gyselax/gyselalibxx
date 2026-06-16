

# Class polar\_foot\_finder\_details::ElementwisePhysicalAdvPhysicalFootFinder

**template &lt;class [**GridR**](structGridR.md), class [**GridTheta**](structGridTheta.md), class IdxRangeOperator, class RThetaAdvectionEvaluator, class AdvecCoefField, class TimeStepper, class LogicalToPhysicalMapping&gt;**



[**ClassList**](annotated.md) **>** [**polar\_foot\_finder\_details**](namespacepolar__foot__finder__details.md) **>** [**ElementwisePhysicalAdvPhysicalFootFinder**](classpolar__foot__finder__details_1_1ElementwisePhysicalAdvPhysicalFootFinder.md)



_GPU-callable functor that finds characteristic feet for physical-space advection with foot-finding in physical_ \((x, y)\) _space._[More...](#detailed-description)

* `#include <physical_advection_physical_foot_finder.hpp>`





































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**ElementwisePhysicalAdvPhysicalFootFinder**](#function-elementwisephysicaladvphysicalfootfinder) (RThetaAdvectionEvaluator const & evaluator\_advection\_field, TimeStepper const & time\_stepper, LogicalToPhysicalMapping logical\_to\_physical, AdvecCoefField const & advection\_field\_coefs, IdxRange&lt; [**GridTheta**](structGridTheta.md) &gt; idx\_range\_theta, Coord&lt; AdvDim1, AdvDim2 &gt; coord\_centre, double dt) <br>_Construct an_ [_**ElementwisePhysicalAdvPhysicalFootFinder**_](classpolar__foot__finder__details_1_1ElementwisePhysicalAdvPhysicalFootFinder.md) _._ |
|  KOKKOS\_FUNCTION CoordRTheta | [**operator()**](#function-operator) (IdxOperator const idx) const<br>_Find the foot of the characteristic at a single grid index._  |




























## Detailed Description


The advection field and foot-finding both operate in the physical Cartesian domain. Constructed by [**ElementwisePhysicalAdvPhysicalFootFinderMem**](classpolar__foot__finder__details_1_1ElementwisePhysicalAdvPhysicalFootFinderMem.md). 


    
## Public Functions Documentation




### function ElementwisePhysicalAdvPhysicalFootFinder 

_Construct an_ [_**ElementwisePhysicalAdvPhysicalFootFinder**_](classpolar__foot__finder__details_1_1ElementwisePhysicalAdvPhysicalFootFinder.md) _._
```C++
inline polar_foot_finder_details::ElementwisePhysicalAdvPhysicalFootFinder::ElementwisePhysicalAdvPhysicalFootFinder (
    RThetaAdvectionEvaluator const & evaluator_advection_field,
    TimeStepper const & time_stepper,
    LogicalToPhysicalMapping logical_to_physical,
    AdvecCoefField const & advection_field_coefs,
    IdxRange< GridTheta > idx_range_theta,
    Coord< AdvDim1, AdvDim2 > coord_centre,
    double dt
) 
```





**Parameters:**


* `evaluator_advection_field` The evaluator for the spline representation of the advection field. 
* `time_stepper` The time integration method. 
* `logical_to_physical` The mapping from the logical domain to the physical domain. 
* `advection_field_coefs` A view of the spline coefficients of the advection field. 
* `idx_range_theta` The poloidal index range, used to wrap the angular coordinate into the periodic domain after each time step. 
* `coord_centre` The coordinate of the polar centre in the physical domain, used to handle the degenerate point at \(r = 0\). 
* `dt` The time step for the characteristic equation. 




        

<hr>



### function operator() 

_Find the foot of the characteristic at a single grid index._ 
```C++
inline KOKKOS_FUNCTION CoordRTheta polar_foot_finder_details::ElementwisePhysicalAdvPhysicalFootFinder::operator() (
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
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/advection/polar_foot_finders/physical_advection_physical_foot_finder.hpp`

