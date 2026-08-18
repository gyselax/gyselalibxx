

# Class polar\_foot\_finder\_details::ElementwisePhysicalAdvPseudoPhysicalFootFinder

**template &lt;class [**GridR**](structGridR.md), class [**GridTheta**](structGridTheta.md), class IdxRangeOperator, class RThetaAdvectionEvaluator, class PseudoPhysicalToAdvectionMapping, class PseudoPhysicalToLogicalMapping, class LogicalToPseudoPhysicalMapping, class AdvecCoefField, class TimeStepper&gt;**



[**ClassList**](annotated.md) **>** [**polar\_foot\_finder\_details**](namespacepolar__foot__finder__details.md) **>** [**ElementwisePhysicalAdvPseudoPhysicalFootFinder**](classpolar__foot__finder__details_1_1ElementwisePhysicalAdvPseudoPhysicalFootFinder.md)



_GPU-callable functor that finds characteristic feet for physical-space advection with foot-finding in pseudo-physical_ \((X_{pC}, Y_{pC})\) _space._[More...](#detailed-description)

* `#include <physical_advection_pseudo_physical_foot_finder.hpp>`





































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**ElementwisePhysicalAdvPseudoPhysicalFootFinder**](#function-elementwisephysicaladvpseudophysicalfootfinder) (RThetaAdvectionEvaluator const & evaluator\_advection\_field, PseudoPhysicalToAdvectionMapping const & pseudo\_physical\_to\_advection, PseudoPhysicalToLogicalMapping const & pseudo\_physical\_to\_logical, LogicalToPseudoPhysicalMapping const & logical\_to\_pseudo\_physical, TimeStepper const & time\_stepper, AdvecCoefField const & advection\_field\_coefs, CoordX\_pcY\_pc coord\_centre, IdxRange&lt; [**GridTheta**](structGridTheta.md) &gt; idx\_range\_theta, double dt) <br>_Construct an_ [_**ElementwisePhysicalAdvPseudoPhysicalFootFinder**_](classpolar__foot__finder__details_1_1ElementwisePhysicalAdvPseudoPhysicalFootFinder.md) _._ |
|  KOKKOS\_FUNCTION CoordRTheta | [**operator()**](#function-operator) (IdxOperator const idx) const<br>_Find the foot of the characteristic at a single grid index._  |




























## Detailed Description


The advection field is expressed in physical Cartesian coordinates and converted to pseudo-physical space. Foot-finding integrates the characteristic equation in pseudo-physical space. Constructed by [**ElementwisePhysicalAdvPseudoPhysicalFootFinderMem**](classpolar__foot__finder__details_1_1ElementwisePhysicalAdvPseudoPhysicalFootFinderMem.md). 


    
## Public Functions Documentation




### function ElementwisePhysicalAdvPseudoPhysicalFootFinder 

_Construct an_ [_**ElementwisePhysicalAdvPseudoPhysicalFootFinder**_](classpolar__foot__finder__details_1_1ElementwisePhysicalAdvPseudoPhysicalFootFinder.md) _._
```C++
inline polar_foot_finder_details::ElementwisePhysicalAdvPseudoPhysicalFootFinder::ElementwisePhysicalAdvPseudoPhysicalFootFinder (
    RThetaAdvectionEvaluator const & evaluator_advection_field,
    PseudoPhysicalToAdvectionMapping const & pseudo_physical_to_advection,
    PseudoPhysicalToLogicalMapping const & pseudo_physical_to_logical,
    LogicalToPseudoPhysicalMapping const & logical_to_pseudo_physical,
    TimeStepper const & time_stepper,
    AdvecCoefField const & advection_field_coefs,
    CoordX_pcY_pc coord_centre,
    IdxRange< GridTheta > idx_range_theta,
    double dt
) 
```





**Parameters:**


* `evaluator_advection_field` The evaluator for the spline representation of the advection field. 
* `pseudo_physical_to_advection` The mapping from the pseudo-physical domain to the advection (physical) domain. 
* `pseudo_physical_to_logical` The mapping from the pseudo-physical domain to the logical domain. 
* `logical_to_pseudo_physical` The mapping from the logical domain to the pseudo-physical domain. 
* `time_stepper` The time integration method. 
* `advection_field_coefs` A view of the spline coefficients of the advection field. 
* `coord_centre` The coordinate of the polar centre in the pseudo-physical domain, used to handle the degenerate point at \(r = 0\). 
* `idx_range_theta` The poloidal index range, used to wrap the angular coordinate into the periodic domain after each time step. 
* `dt` The time step for the characteristic equation. 




        

<hr>



### function operator() 

_Find the foot of the characteristic at a single grid index._ 
```C++
inline KOKKOS_FUNCTION CoordRTheta polar_foot_finder_details::ElementwisePhysicalAdvPseudoPhysicalFootFinder::operator() (
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
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/advection/polar_foot_finders/physical_advection_pseudo_physical_foot_finder.hpp`

