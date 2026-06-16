

# Class polar\_foot\_finder\_details::ElementwisePhysicalAdvPhysicalFootFinderMem

**template &lt;class [**GridR**](structGridR.md), class [**GridTheta**](structGridTheta.md), class IdxRangeOperator, class RThetaAdvectionEvaluator, class AdvecCoefFieldMem, class TimeStepperBuilder, class LogicalToPhysicalMapping&gt;**



[**ClassList**](annotated.md) **>** [**polar\_foot\_finder\_details**](namespacepolar__foot__finder__details.md) **>** [**ElementwisePhysicalAdvPhysicalFootFinderMem**](classpolar__foot__finder__details_1_1ElementwisePhysicalAdvPhysicalFootFinderMem.md)



_Owning manager for_ [_**ElementwisePhysicalAdvPhysicalFootFinder**_](classpolar__foot__finder__details_1_1ElementwisePhysicalAdvPhysicalFootFinder.md) _._[More...](#detailed-description)

* `#include <physical_advection_physical_foot_finder.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef [**ElementwisePhysicalAdvPhysicalFootFinder**](classpolar__foot__finder__details_1_1ElementwisePhysicalAdvPhysicalFootFinder.md)&lt; [**GridR**](structGridR.md), [**GridTheta**](structGridTheta.md), IdxRangeOperator, RThetaAdvectionEvaluator, typename AdvecCoefFieldMem::view\_type, TimeStepper, LogicalToPhysicalMapping &gt; | [**GPUCompat**](#typedef-gpucompat)  <br>_The non-owning operator that can be used on GPU._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**ElementwisePhysicalAdvPhysicalFootFinderMem**](#function-elementwisephysicaladvphysicalfootfindermem) (RThetaAdvectionEvaluator const & evaluator\_advection\_field, LogicalToPhysicalMapping const & logical\_to\_physical, TimeStepperBuilder const & time\_stepper\_builder, AdvecCoefFieldMem && advection\_field\_coefs, Coord&lt; AdvDim1, AdvDim2 &gt; coord\_centre, IdxRange&lt; [**GridTheta**](structGridTheta.md) &gt; idx\_range\_theta) <br>_Construct an_ [_**ElementwisePhysicalAdvPhysicalFootFinderMem**_](classpolar__foot__finder__details_1_1ElementwisePhysicalAdvPhysicalFootFinderMem.md) _._ |
|  [**GPUCompat**](classpolar__foot__finder__details_1_1ElementwisePhysicalAdvPhysicalFootFinderMem.md#typedef-gpucompat) | [**operator()**](#function-operator) (double dt) <br>_Create an_ [_**ElementwisePhysicalAdvPhysicalFootFinder**_](classpolar__foot__finder__details_1_1ElementwisePhysicalAdvPhysicalFootFinder.md) _for the given time step._ |




























## Detailed Description


Holds the spline coefficient field and the preallocated time stepper. Call `operator()(dt)` to obtain a GPU-copyable [**ElementwisePhysicalAdvPhysicalFootFinder**](classpolar__foot__finder__details_1_1ElementwisePhysicalAdvPhysicalFootFinder.md) configured for time step `dt`. 


    
## Public Types Documentation




### typedef GPUCompat 

_The non-owning operator that can be used on GPU._ 
```C++
using polar_foot_finder_details::ElementwisePhysicalAdvPhysicalFootFinderMem< GridR, GridTheta, IdxRangeOperator, RThetaAdvectionEvaluator, AdvecCoefFieldMem, TimeStepperBuilder, LogicalToPhysicalMapping >::GPUCompat =  ElementwisePhysicalAdvPhysicalFootFinder< GridR, GridTheta, IdxRangeOperator, RThetaAdvectionEvaluator, typename AdvecCoefFieldMem::view_type, TimeStepper, LogicalToPhysicalMapping>;
```




<hr>
## Public Functions Documentation




### function ElementwisePhysicalAdvPhysicalFootFinderMem 

_Construct an_ [_**ElementwisePhysicalAdvPhysicalFootFinderMem**_](classpolar__foot__finder__details_1_1ElementwisePhysicalAdvPhysicalFootFinderMem.md) _._
```C++
inline polar_foot_finder_details::ElementwisePhysicalAdvPhysicalFootFinderMem::ElementwisePhysicalAdvPhysicalFootFinderMem (
    RThetaAdvectionEvaluator const & evaluator_advection_field,
    LogicalToPhysicalMapping const & logical_to_physical,
    TimeStepperBuilder const & time_stepper_builder,
    AdvecCoefFieldMem && advection_field_coefs,
    Coord< AdvDim1, AdvDim2 > coord_centre,
    IdxRange< GridTheta > idx_range_theta
) 
```





**Parameters:**


* `evaluator_advection_field` The evaluator for the spline representation of the advection field. 
* `logical_to_physical` The mapping from the logical domain to the physical domain. 
* `time_stepper_builder` The factory used to preallocate the time integration method. 
* `advection_field_coefs` The spline coefficients of the advection field. Ownership is transferred in. 
* `coord_centre` The coordinate of the polar centre in the physical domain, used to handle the degenerate point at \(r = 0\). 
* `idx_range_theta` The poloidal index range, used to wrap the angular coordinate into the periodic domain after each time step. 




        

<hr>



### function operator() 

_Create an_ [_**ElementwisePhysicalAdvPhysicalFootFinder**_](classpolar__foot__finder__details_1_1ElementwisePhysicalAdvPhysicalFootFinder.md) _for the given time step._
```C++
inline GPUCompat polar_foot_finder_details::ElementwisePhysicalAdvPhysicalFootFinderMem::operator() (
    double dt
) 
```



Returns a non-owning [**ElementwisePhysicalAdvPhysicalFootFinder**](classpolar__foot__finder__details_1_1ElementwisePhysicalAdvPhysicalFootFinder.md) that holds views of the stored spline coefficients and is configured for time step \(dt\). The returned object can be copied to and called from the device.




**Parameters:**


* `dt` The time step for the characteristic equation.



**Returns:**

A view-based [**ElementwisePhysicalAdvPhysicalFootFinder**](classpolar__foot__finder__details_1_1ElementwisePhysicalAdvPhysicalFootFinder.md) for the given time step. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/advection/polar_foot_finders/physical_advection_physical_foot_finder.hpp`

