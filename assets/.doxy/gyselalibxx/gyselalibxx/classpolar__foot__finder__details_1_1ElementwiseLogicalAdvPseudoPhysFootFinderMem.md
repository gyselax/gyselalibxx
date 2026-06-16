

# Class polar\_foot\_finder\_details::ElementwiseLogicalAdvPseudoPhysFootFinderMem

**template &lt;class [**GridR**](structGridR.md), class [**GridTheta**](structGridTheta.md), class IdxRangeOperator, class RThetaAdvectionEvaluator, class AdvecCoefFieldMem, class TimeStepperBuilder, class LogicalToPseudoPhysicalMapping&gt;**



[**ClassList**](annotated.md) **>** [**polar\_foot\_finder\_details**](namespacepolar__foot__finder__details.md) **>** [**ElementwiseLogicalAdvPseudoPhysFootFinderMem**](classpolar__foot__finder__details_1_1ElementwiseLogicalAdvPseudoPhysFootFinderMem.md)



_Owning manager for_ [_**ElementwiseLogicalAdvPseudoPhysFootFinder**_](classpolar__foot__finder__details_1_1ElementwiseLogicalAdvPseudoPhysFootFinder.md) _._[More...](#detailed-description)

* `#include <logical_advection_pseudo_physical_foot_finder.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef [**ElementwiseLogicalAdvPseudoPhysFootFinder**](classpolar__foot__finder__details_1_1ElementwiseLogicalAdvPseudoPhysFootFinder.md)&lt; [**GridR**](structGridR.md), [**GridTheta**](structGridTheta.md), IdxRangeOperator, RThetaAdvectionEvaluator, typename AdvecCoefFieldMem::view\_type, TimeStepper, LogicalToPseudoPhysicalMapping &gt; | [**GPUCompat**](#typedef-gpucompat)  <br>_The non-owning operator that can be used on GPU._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**ElementwiseLogicalAdvPseudoPhysFootFinderMem**](#function-elementwiselogicaladvpseudophysfootfindermem-12) (RThetaAdvectionEvaluator const & evaluator\_advection\_field, LogicalToPhysicalMapping const & logical\_to\_physical, TimeStepperBuilder const & time\_stepper\_builder, AdvecCoefFieldMem && advection\_field\_coefs, Coord&lt; X\_pc, Y\_pc &gt; coord\_centre, IdxRange&lt; [**GridTheta**](structGridTheta.md) &gt; idx\_range\_theta) <br>_Construct an_ [_**ElementwiseLogicalAdvPseudoPhysFootFinderMem**_](classpolar__foot__finder__details_1_1ElementwiseLogicalAdvPseudoPhysFootFinderMem.md) _._ |
|   | [**ElementwiseLogicalAdvPseudoPhysFootFinderMem**](#function-elementwiselogicaladvpseudophysfootfindermem-12) (RThetaAdvectionEvaluator const & evaluator\_advection\_field, LogicalToPhysicalMapping const & logical\_to\_physical, TimeStepperBuilder const & time\_stepper\_builder, AdvecCoefFieldMem && advection\_field\_coefs, Coord&lt; X\_pc, Y\_pc &gt; coord\_centre, IdxRange&lt; [**GridTheta**](structGridTheta.md) &gt; idx\_range\_theta) <br>_Construct an_ [_**ElementwiseLogicalAdvPseudoPhysFootFinderMem**_](classpolar__foot__finder__details_1_1ElementwiseLogicalAdvPseudoPhysFootFinderMem.md) _when the logical-to-physical mapping does not share the pseudo-physical Cartesian tags._ |
|  [**GPUCompat**](classpolar__foot__finder__details_1_1ElementwiseLogicalAdvPseudoPhysFootFinderMem.md#typedef-gpucompat) | [**operator()**](#function-operator) (double dt) <br>_Create an_ [_**ElementwiseLogicalAdvPseudoPhysFootFinder**_](classpolar__foot__finder__details_1_1ElementwiseLogicalAdvPseudoPhysFootFinder.md) _for the given time step._ |




























## Detailed Description


Holds the spline coefficient field and the preallocated time stepper. Call `operator()(dt)` to obtain a GPU-copyable [**ElementwiseLogicalAdvPseudoPhysFootFinder**](classpolar__foot__finder__details_1_1ElementwiseLogicalAdvPseudoPhysFootFinder.md) configured for time step `dt`. 


    
## Public Types Documentation




### typedef GPUCompat 

_The non-owning operator that can be used on GPU._ 
```C++
using polar_foot_finder_details::ElementwiseLogicalAdvPseudoPhysFootFinderMem< GridR, GridTheta, IdxRangeOperator, RThetaAdvectionEvaluator, AdvecCoefFieldMem, TimeStepperBuilder, LogicalToPseudoPhysicalMapping >::GPUCompat =  ElementwiseLogicalAdvPseudoPhysFootFinder< GridR, GridTheta, IdxRangeOperator, RThetaAdvectionEvaluator, typename AdvecCoefFieldMem::view_type, TimeStepper, LogicalToPseudoPhysicalMapping>;
```




<hr>
## Public Functions Documentation




### function ElementwiseLogicalAdvPseudoPhysFootFinderMem [1/2]

_Construct an_ [_**ElementwiseLogicalAdvPseudoPhysFootFinderMem**_](classpolar__foot__finder__details_1_1ElementwiseLogicalAdvPseudoPhysFootFinderMem.md) _._
```C++
template<class LogicalToPhysicalMapping, std::enable_if_t<(std::is_same_v< X_pc, typename LogicalToPhysicalMapping::cartesian_tag_x >)&&(std::is_same_v< Y_pc, typename LogicalToPhysicalMapping::cartesian_tag_y >), bool >>
inline polar_foot_finder_details::ElementwiseLogicalAdvPseudoPhysFootFinderMem::ElementwiseLogicalAdvPseudoPhysFootFinderMem (
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
* `logical_to_physical` The mapping from the logical domain to the physical domain, used to derive the pseudo-physical mappings. 
* `time_stepper_builder` The factory used to preallocate the time integration method. 
* `advection_field_coefs` The spline coefficients of the advection field. Ownership is transferred in. 
* `coord_centre` The coordinate of the polar centre in the pseudo-physical domain, used to handle the degenerate point at \(r = 0\). 
* `idx_range_theta` The poloidal index range, used to wrap the angular coordinate into the periodic domain after each time step. 




        

<hr>



### function ElementwiseLogicalAdvPseudoPhysFootFinderMem [1/2]

_Construct an_ [_**ElementwiseLogicalAdvPseudoPhysFootFinderMem**_](classpolar__foot__finder__details_1_1ElementwiseLogicalAdvPseudoPhysFootFinderMem.md) _when the logical-to-physical mapping does not share the pseudo-physical Cartesian tags._
```C++
template<class LogicalToPhysicalMapping, std::enable_if_t< !((std::is_same_v< X_pc, typename LogicalToPhysicalMapping::cartesian_tag_x >)&&(std::is_same_v< Y_pc, typename LogicalToPhysicalMapping::cartesian_tag_y >)), bool >>
inline polar_foot_finder_details::ElementwiseLogicalAdvPseudoPhysFootFinderMem::ElementwiseLogicalAdvPseudoPhysFootFinderMem (
    RThetaAdvectionEvaluator const & evaluator_advection_field,
    LogicalToPhysicalMapping const & logical_to_physical,
    TimeStepperBuilder const & time_stepper_builder,
    AdvecCoefFieldMem && advection_field_coefs,
    Coord< X_pc, Y_pc > coord_centre,
    IdxRange< GridTheta > idx_range_theta
) 
```



In this case `logical_to_physical` is ignored and the pseudo-physical mapping is constructed directly from `coord_centre` using a default circular mapping.




**Parameters:**


* `evaluator_advection_field` The evaluator for the spline representation of the advection field. 
* `logical_to_physical` Unused; accepted for interface uniformity with the primary constructor. 
* `time_stepper_builder` The factory used to preallocate the time integration method. 
* `advection_field_coefs` The spline coefficients of the advection field. Ownership is transferred in. 
* `coord_centre` The coordinate of the polar centre in the pseudo-physical domain, used to initialise the pseudo-physical mappings and handle the degenerate point at \(r = 0\). 
* `idx_range_theta` The poloidal index range, used to wrap the angular coordinate into the periodic domain after each time step. 




        

<hr>



### function operator() 

_Create an_ [_**ElementwiseLogicalAdvPseudoPhysFootFinder**_](classpolar__foot__finder__details_1_1ElementwiseLogicalAdvPseudoPhysFootFinder.md) _for the given time step._
```C++
inline GPUCompat polar_foot_finder_details::ElementwiseLogicalAdvPseudoPhysFootFinderMem::operator() (
    double dt
) 
```



Returns a non-owning [**ElementwiseLogicalAdvPseudoPhysFootFinder**](classpolar__foot__finder__details_1_1ElementwiseLogicalAdvPseudoPhysFootFinder.md) that holds views of the stored spline coefficients and is configured for time step \(dt\). The returned object can be copied to and called from the device.




**Parameters:**


* `dt` The time step for the characteristic equation.



**Returns:**

A view-based [**ElementwiseLogicalAdvPseudoPhysFootFinder**](classpolar__foot__finder__details_1_1ElementwiseLogicalAdvPseudoPhysFootFinder.md) for the given time step. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/advection/polar_foot_finders/logical_advection_pseudo_physical_foot_finder.hpp`

