

# Class ElementwiseSplinePolarFootFinder

**template &lt;class [**GridR**](structGridR.md), class [**GridTheta**](structGridTheta.md), class IdxRangeOperator, class SplineRThetaEvaluatorAdvection, class PseudoPhysicalToAdvectionMapping, class PseudoPhysicalToLogicalMapping, class LogicalToPseudoPhysicalMapping, class AdvecCoefField, class TimeStepper&gt;**



[**ClassList**](annotated.md) **>** [**ElementwiseSplinePolarFootFinder**](classElementwiseSplinePolarFootFinder.md)



_A device-callable functor that finds the foot of the characteristic at a single grid point in polar coordinates using spline interpolation._ [More...](#detailed-description)

* `#include <spline_polar_foot_finder.hpp>`





































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**ElementwiseSplinePolarFootFinder**](#function-elementwisesplinepolarfootfinder) (SplineRThetaEvaluatorAdvection const & evaluator\_advection\_field, PseudoPhysicalToAdvectionMapping const & pseudo\_physical\_to\_advection, PseudoPhysicalToLogicalMapping const & pseudo\_physical\_to\_logical, LogicalToPseudoPhysicalMapping const & logical\_to\_pseudo\_physical, TimeStepper const & time\_stepper, AdvecCoefField const & advection\_field\_coefs, CoordX\_pcY\_pc coord\_centre, IdxRange&lt; [**GridTheta**](structGridTheta.md) &gt; idx\_range\_theta, double dt) <br>_Construct an_ [_**ElementwiseSplinePolarFootFinder**_](classElementwiseSplinePolarFootFinder.md) _._ |
|  KOKKOS\_FUNCTION CoordRTheta | [**operator()**](#function-operator) (IdxOperator const idx) const<br>_Find the foot of the characteristic at a single grid index._  |




























## Detailed Description


This class holds non-owning views of all the objects needed to solve the characteristic equation at a single index, and so it can be copied to and called from the device. It is the innermost computation unit of [**SplinePolarFootFinder**](classSplinePolarFootFinder.md). The calculation of the foot is carried out in Cartesian or pseudo-Cartesian coordinates to avoir problems related to the O-point, namely:
* The advection field is undefined in polar coordinates at the O-point.
* Calculating the foot of a characteristic crossing the O-point could lead to negative radial values.






**Template parameters:**


* [**GridR**](structGridR.md) The discrete radial dimension. 
* [**GridTheta**](structGridTheta.md) The discrete poloidal dimension. 
* `IdxRangeOperator` The full index range over which the operator acts (may include batch dimensions). 
* `SplineRThetaEvaluatorAdvection` The evaluator used to evaluate the spline representation of the advection field. 
* `PseudoPhysicalToAdvectionMapping` A mapping from the pseudo-physical domain to the domain where the advection field is defined. 
* `PseudoPhysicalToLogicalMapping` A mapping from the pseudo-physical domain to the logical domain. 
* `LogicalToPseudoPhysicalMapping` A mapping from the logical domain to the pseudo-physical domain. 
* `AdvecCoefField` A non-owning field (view) type holding the spline coefficients of the advection field. 
* `TimeStepper` The time integration method used to solve the characteristic equation.



**See also:** [**ElementwiseSplinePolarFootFinderMem**](classElementwiseSplinePolarFootFinderMem.md) 


**See also:** [**SplinePolarFootFinder**](classSplinePolarFootFinder.md) 



    
## Public Functions Documentation




### function ElementwiseSplinePolarFootFinder 

_Construct an_ [_**ElementwiseSplinePolarFootFinder**_](classElementwiseSplinePolarFootFinder.md) _._
```C++
inline ElementwiseSplinePolarFootFinder::ElementwiseSplinePolarFootFinder (
    SplineRThetaEvaluatorAdvection const & evaluator_advection_field,
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
* `pseudo_physical_to_advection` The mapping from the pseudo-physical domain to the advection field domain. 
* `pseudo_physical_to_logical` The mapping from the pseudo-physical domain to the logical domain. 
* `logical_to_pseudo_physical` The mapping from the logical domain to the pseudo-physical domain. 
* `time_stepper` The time integration method. 
* `advection_field_coefs` A non-owning view of the pre-built spline coefficients of the advection field. 
* `coord_centre` The coordinate of the polar centre in the pseudo-physical domain, used to handle the degenerate point at \(r = 0\). 
* `idx_range_theta` The poloidal index range, used to wrap the angular coordinate into the periodic domain after each time step. 
* `dt` The time step for the characteristic equation. 




        

<hr>



### function operator() 

_Find the foot of the characteristic at a single grid index._ 
```C++
inline KOKKOS_FUNCTION CoordRTheta ElementwiseSplinePolarFootFinder::operator() (
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
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/advection/spline_polar_foot_finder.hpp`

