

# Class ElementwiseSplinePolarFootFinderMem

**template &lt;class [**GridR**](structGridR.md), class [**GridTheta**](structGridTheta.md), class IdxRangeOperator, class SplineRThetaEvaluatorAdvection, class PseudoPhysicalToAdvectionMapping, class PseudoPhysicalToLogicalMapping, class LogicalToPseudoPhysicalMapping, class AdvecCoefFieldMem, class TimeStepper&gt;**



[**ClassList**](annotated.md) **>** [**ElementwiseSplinePolarFootFinderMem**](classElementwiseSplinePolarFootFinderMem.md)



_The owning counterpart to_ [_**ElementwiseSplinePolarFootFinder**_](classElementwiseSplinePolarFootFinder.md) _._[More...](#detailed-description)

* `#include <spline_polar_foot_finder.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef [**ElementwiseSplinePolarFootFinder**](classElementwiseSplinePolarFootFinder.md)&lt; [**GridR**](structGridR.md), [**GridTheta**](structGridTheta.md), IdxRangeOperator, SplineRThetaEvaluatorAdvection, PseudoPhysicalToAdvectionMapping, PseudoPhysicalToLogicalMapping, LogicalToPseudoPhysicalMapping, typename AdvecCoefFieldMem::view\_type, TimeStepper &gt; | [**GPUCompat**](#typedef-gpucompat)  <br>_The non-owning operator that can be used on GPU._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**ElementwiseSplinePolarFootFinderMem**](#function-elementwisesplinepolarfootfindermem) (SplineRThetaEvaluatorAdvection const & evaluator\_advection\_field, PseudoPhysicalToAdvectionMapping const & pseudo\_physical\_to\_advection, PseudoPhysicalToLogicalMapping const & pseudo\_physical\_to\_logical, LogicalToPseudoPhysicalMapping const & logical\_to\_pseudo\_physical, TimeStepper const & time\_stepper, AdvecCoefFieldMem && advection\_field\_coefs, Coord&lt; X\_pc, Y\_pc &gt; coord\_centre, IdxRange&lt; [**GridTheta**](structGridTheta.md) &gt; idx\_range\_theta) <br>_Construct an_ [_**ElementwiseSplinePolarFootFinderMem**_](classElementwiseSplinePolarFootFinderMem.md) _._ |
|  [**GPUCompat**](classElementwiseSplinePolarFootFinderMem.md#typedef-gpucompat) | [**operator()**](#function-operator) (double dt) <br>_Create an_ [_**ElementwiseSplinePolarFootFinder**_](classElementwiseSplinePolarFootFinder.md) _for the given time step._ |




























## Detailed Description


Allocates and stores the spline coefficients of the advection field on the appropriate memory space. Calling `operator()(dt)` returns a non-owning [**ElementwiseSplinePolarFootFinder**](classElementwiseSplinePolarFootFinder.md) configured for the given time step, which can then be copied to and called from the device.




**Template parameters:**


* [**GridR**](structGridR.md) The discrete radial dimension. 
* [**GridTheta**](structGridTheta.md) The discrete poloidal dimension. 
* `IdxRangeOperator` The full index range over which the operator acts (may include batch dimensions). 
* `SplineRThetaEvaluatorAdvection` The evaluator used to evaluate the spline representation of the advection field. 
* `PseudoPhysicalToAdvectionMapping` A mapping from the pseudo-physical domain to the domain where the advection field is defined. 
* `PseudoPhysicalToLogicalMapping` A mapping from the pseudo-physical domain to the logical domain. 
* `LogicalToPseudoPhysicalMapping` A mapping from the logical domain to the pseudo-physical domain. 
* `AdvecCoefFieldMem` An owning field type holding the spline coefficients of the advection field. 
* `TimeStepper` The time integration method used to solve the characteristic equation.



**See also:** [**ElementwiseSplinePolarFootFinder**](classElementwiseSplinePolarFootFinder.md) 


**See also:** [**SplinePolarFootFinder**](classSplinePolarFootFinder.md) 



    
## Public Types Documentation




### typedef GPUCompat 

_The non-owning operator that can be used on GPU._ 
```C++
using ElementwiseSplinePolarFootFinderMem< GridR, GridTheta, IdxRangeOperator, SplineRThetaEvaluatorAdvection, PseudoPhysicalToAdvectionMapping, PseudoPhysicalToLogicalMapping, LogicalToPseudoPhysicalMapping, AdvecCoefFieldMem, TimeStepper >::GPUCompat =  ElementwiseSplinePolarFootFinder< GridR, GridTheta, IdxRangeOperator, SplineRThetaEvaluatorAdvection, PseudoPhysicalToAdvectionMapping, PseudoPhysicalToLogicalMapping, LogicalToPseudoPhysicalMapping, typename AdvecCoefFieldMem::view_type, TimeStepper>;
```




<hr>
## Public Functions Documentation




### function ElementwiseSplinePolarFootFinderMem 

_Construct an_ [_**ElementwiseSplinePolarFootFinderMem**_](classElementwiseSplinePolarFootFinderMem.md) _._
```C++
inline ElementwiseSplinePolarFootFinderMem::ElementwiseSplinePolarFootFinderMem (
    SplineRThetaEvaluatorAdvection const & evaluator_advection_field,
    PseudoPhysicalToAdvectionMapping const & pseudo_physical_to_advection,
    PseudoPhysicalToLogicalMapping const & pseudo_physical_to_logical,
    LogicalToPseudoPhysicalMapping const & logical_to_pseudo_physical,
    TimeStepper const & time_stepper,
    AdvecCoefFieldMem && advection_field_coefs,
    Coord< X_pc, Y_pc > coord_centre,
    IdxRange< GridTheta > idx_range_theta
) 
```





**Parameters:**


* `evaluator_advection_field` The evaluator for the spline representation of the advection field. 
* `pseudo_physical_to_advection` The mapping from the pseudo-physical domain to the advection field domain. 
* `pseudo_physical_to_logical` The mapping from the pseudo-physical domain to the logical domain. 
* `logical_to_pseudo_physical` The mapping from the logical domain to the pseudo-physical domain. 
* `time_stepper` The time integration method. 
* `advection_field_coefs` The spline coefficients of the advection field. Ownership is transferred in. 
* `coord_centre` The coordinate of the polar centre in the pseudo-physical domain, used to handle the degenerate point at \(r = 0\). 
* `idx_range_theta` The poloidal index range, used to wrap the angular coordinate into the periodic domain after each time step. 




        

<hr>



### function operator() 

_Create an_ [_**ElementwiseSplinePolarFootFinder**_](classElementwiseSplinePolarFootFinder.md) _for the given time step._
```C++
inline GPUCompat ElementwiseSplinePolarFootFinderMem::operator() (
    double dt
) 
```



Returns a non-owning [**ElementwiseSplinePolarFootFinder**](classElementwiseSplinePolarFootFinder.md) that holds views of the stored spline coefficients and is configured for time step \(dt\). The returned object can be copied to and called from the device.




**Parameters:**


* `dt` The time step for the characteristic equation.



**Returns:**

A view-based [**ElementwiseSplinePolarFootFinder**](classElementwiseSplinePolarFootFinder.md) for the given time step. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/advection/spline_polar_foot_finder.hpp`

