

# Class polar\_foot\_finder\_details::ElementwiseLogicalAdvLogicalFootFinderMem

**template &lt;class [**GridR**](structGridR.md), class [**GridTheta**](structGridTheta.md), class IdxRangeOperator, class RThetaAdvectionEvaluator, class AdvecCoefFieldMem, class TimeStepperBuilder&gt;**



[**ClassList**](annotated.md) **>** [**polar\_foot\_finder\_details**](namespacepolar__foot__finder__details.md) **>** [**ElementwiseLogicalAdvLogicalFootFinderMem**](classpolar__foot__finder__details_1_1ElementwiseLogicalAdvLogicalFootFinderMem.md)



_Owning manager for_ [_**ElementwiseLogicalAdvLogicalFootFinder**_](classpolar__foot__finder__details_1_1ElementwiseLogicalAdvLogicalFootFinder.md) _._[More...](#detailed-description)

* `#include <logical_advection_logical_foot_finder.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef [**ElementwiseLogicalAdvLogicalFootFinder**](classpolar__foot__finder__details_1_1ElementwiseLogicalAdvLogicalFootFinder.md)&lt; [**GridR**](structGridR.md), [**GridTheta**](structGridTheta.md), IdxRangeOperator, RThetaAdvectionEvaluator, typename AdvecCoefFieldMem::view\_type, TimeStepper &gt; | [**GPUCompat**](#typedef-gpucompat)  <br>_The non-owning operator that can be used on GPU._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**ElementwiseLogicalAdvLogicalFootFinderMem**](#function-elementwiselogicaladvlogicalfootfindermem) (RThetaAdvectionEvaluator const & evaluator\_advection\_field, LogicalToPhysicalMapping const & logical\_to\_physical, TimeStepperBuilder const & time\_stepper\_builder, AdvecCoefFieldMem && advection\_field\_coefs, Coord&lt; X\_pc, Y\_pc &gt; coord\_centre, IdxRange&lt; [**GridTheta**](structGridTheta.md) &gt; idx\_range\_theta) <br>_Construct an_ [_**ElementwiseLogicalAdvLogicalFootFinderMem**_](classpolar__foot__finder__details_1_1ElementwiseLogicalAdvLogicalFootFinderMem.md) _._ |
|  [**GPUCompat**](classpolar__foot__finder__details_1_1ElementwiseLogicalAdvLogicalFootFinderMem.md#typedef-gpucompat) | [**operator()**](#function-operator) (double dt) <br>_Create an_ [_**ElementwiseLogicalAdvLogicalFootFinder**_](classpolar__foot__finder__details_1_1ElementwiseLogicalAdvLogicalFootFinder.md) _for the given time step._ |




























## Detailed Description


Holds the spline coefficient field and the preallocated time stepper. Call `operator()(dt)` to obtain a GPU-copyable [**ElementwiseLogicalAdvLogicalFootFinder**](classpolar__foot__finder__details_1_1ElementwiseLogicalAdvLogicalFootFinder.md) configured for time step `dt`. 


    
## Public Types Documentation




### typedef GPUCompat 

_The non-owning operator that can be used on GPU._ 
```C++
using polar_foot_finder_details::ElementwiseLogicalAdvLogicalFootFinderMem< GridR, GridTheta, IdxRangeOperator, RThetaAdvectionEvaluator, AdvecCoefFieldMem, TimeStepperBuilder >::GPUCompat =  ElementwiseLogicalAdvLogicalFootFinder< GridR, GridTheta, IdxRangeOperator, RThetaAdvectionEvaluator, typename AdvecCoefFieldMem::view_type, TimeStepper>;
```




<hr>
## Public Functions Documentation




### function ElementwiseLogicalAdvLogicalFootFinderMem 

_Construct an_ [_**ElementwiseLogicalAdvLogicalFootFinderMem**_](classpolar__foot__finder__details_1_1ElementwiseLogicalAdvLogicalFootFinderMem.md) _._
```C++
template<class LogicalToPhysicalMapping, class X_pc, class Y_pc>
inline polar_foot_finder_details::ElementwiseLogicalAdvLogicalFootFinderMem::ElementwiseLogicalAdvLogicalFootFinderMem (
    RThetaAdvectionEvaluator const & evaluator_advection_field,
    LogicalToPhysicalMapping const & logical_to_physical,
    TimeStepperBuilder const & time_stepper_builder,
    AdvecCoefFieldMem && advection_field_coefs,
    Coord< X_pc, Y_pc > coord_centre,
    IdxRange< GridTheta > idx_range_theta
) 
```





**Parameters:**


* `evaluator_advection_field` The evaluator for the spline representation of the advection field. 
* `logical_to_physical` Unused for the logical/logical case; accepted for interface uniformity with other [**ElementwiseChoice**](structpolar__foot__finder__details_1_1ElementwiseChoice.md) specialisations. 
* `time_stepper_builder` The factory used to preallocate the time integration method. 
* `advection_field_coefs` The spline coefficients of the advection field. Ownership is transferred in. 
* `coord_centre` Unused for the logical/logical case; accepted for interface uniformity with other [**ElementwiseChoice**](structpolar__foot__finder__details_1_1ElementwiseChoice.md) specialisations. 
* `idx_range_theta` The poloidal index range, used to wrap the angular coordinate into the periodic domain after each time step. 




        

<hr>



### function operator() 

_Create an_ [_**ElementwiseLogicalAdvLogicalFootFinder**_](classpolar__foot__finder__details_1_1ElementwiseLogicalAdvLogicalFootFinder.md) _for the given time step._
```C++
inline GPUCompat polar_foot_finder_details::ElementwiseLogicalAdvLogicalFootFinderMem::operator() (
    double dt
) 
```



Returns a non-owning [**ElementwiseLogicalAdvLogicalFootFinder**](classpolar__foot__finder__details_1_1ElementwiseLogicalAdvLogicalFootFinder.md) that holds views of the stored spline coefficients and is configured for time step \(dt\). The returned object can be copied to and called from the device.




**Parameters:**


* `dt` The time step for the characteristic equation.



**Returns:**

A view-based [**ElementwiseLogicalAdvLogicalFootFinder**](classpolar__foot__finder__details_1_1ElementwiseLogicalAdvLogicalFootFinder.md) for the given time step. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/advection/polar_foot_finders/logical_advection_logical_foot_finder.hpp`

