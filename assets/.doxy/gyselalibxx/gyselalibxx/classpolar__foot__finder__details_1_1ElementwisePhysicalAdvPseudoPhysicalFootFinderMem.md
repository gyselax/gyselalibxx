

# Class polar\_foot\_finder\_details::ElementwisePhysicalAdvPseudoPhysicalFootFinderMem

**template &lt;class [**GridR**](structGridR.md), class [**GridTheta**](structGridTheta.md), class IdxRangeOperator, class RThetaAdvectionEvaluator, class AdvecCoefFieldMem, class TimeStepperBuilder, concepts::Mapping LogicalToPhysicalMapping&gt;**



[**ClassList**](annotated.md) **>** [**polar\_foot\_finder\_details**](namespacepolar__foot__finder__details.md) **>** [**ElementwisePhysicalAdvPseudoPhysicalFootFinderMem**](classpolar__foot__finder__details_1_1ElementwisePhysicalAdvPseudoPhysicalFootFinderMem.md)



_Owning manager for_ [_**ElementwisePhysicalAdvPseudoPhysicalFootFinder**_](classpolar__foot__finder__details_1_1ElementwisePhysicalAdvPseudoPhysicalFootFinder.md) _._[More...](#detailed-description)

* `#include <physical_advection_pseudo_physical_foot_finder.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef [**ElementwisePhysicalAdvPseudoPhysicalFootFinder**](classpolar__foot__finder__details_1_1ElementwisePhysicalAdvPseudoPhysicalFootFinder.md)&lt; [**GridR**](structGridR.md), [**GridTheta**](structGridTheta.md), IdxRangeOperator, RThetaAdvectionEvaluator, [**PseudoPhysicalToAdvectionMapping**](classCombinedMapping.md), [**PseudoPhysicalToLogicalMapping**](classCartesianToCircular.md), [**LogicalToPseudoPhysicalMapping**](classCircularToCartesian.md), typename AdvecCoefFieldMem::view\_type, TimeStepper &gt; | [**GPUCompat**](#typedef-gpucompat)  <br>_The non-owning operator that can be used on GPU._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**ElementwisePhysicalAdvPseudoPhysicalFootFinderMem**](#function-elementwisephysicaladvpseudophysicalfootfindermem) (RThetaAdvectionEvaluator const & evaluator\_advection\_field, LogicalToPhysicalMapping const & logical\_to\_physical, TimeStepperBuilder const & time\_stepper\_builder, AdvecCoefFieldMem && advection\_field\_coefs, Coord&lt; [**X\_pC**](structX__pC.md), [**Y\_pC**](structY__pC.md) &gt; coord\_centre, IdxRange&lt; [**GridTheta**](structGridTheta.md) &gt; idx\_range\_theta, double epsilon=1e-12) <br>_Construct an_ [_**ElementwisePhysicalAdvPseudoPhysicalFootFinderMem**_](classpolar__foot__finder__details_1_1ElementwisePhysicalAdvPseudoPhysicalFootFinderMem.md) _._ |
|  [**GPUCompat**](classpolar__foot__finder__details_1_1ElementwisePhysicalAdvPseudoPhysicalFootFinderMem.md#typedef-gpucompat) | [**operator()**](#function-operator) (double dt) <br>_Create an_ [_**ElementwisePhysicalAdvPseudoPhysicalFootFinder**_](classpolar__foot__finder__details_1_1ElementwisePhysicalAdvPseudoPhysicalFootFinder.md) _for the given time step._ |




























## Detailed Description


Builds the `CombinedMapping` (physical-to-advection), holds the spline coefficient field, and preallocates the time stepper. Call `operator()(dt)` to obtain a GPU-copyable [**ElementwisePhysicalAdvPseudoPhysicalFootFinder**](classpolar__foot__finder__details_1_1ElementwisePhysicalAdvPseudoPhysicalFootFinder.md) configured for time step `dt`. 


    
## Public Types Documentation




### typedef GPUCompat 

_The non-owning operator that can be used on GPU._ 
```C++
using polar_foot_finder_details::ElementwisePhysicalAdvPseudoPhysicalFootFinderMem< GridR, GridTheta, IdxRangeOperator, RThetaAdvectionEvaluator, AdvecCoefFieldMem, TimeStepperBuilder, LogicalToPhysicalMapping >::GPUCompat =  ElementwisePhysicalAdvPseudoPhysicalFootFinder< GridR, GridTheta, IdxRangeOperator, RThetaAdvectionEvaluator, PseudoPhysicalToAdvectionMapping, PseudoPhysicalToLogicalMapping, LogicalToPseudoPhysicalMapping, typename AdvecCoefFieldMem::view_type, TimeStepper>;
```




<hr>
## Public Functions Documentation




### function ElementwisePhysicalAdvPseudoPhysicalFootFinderMem 

_Construct an_ [_**ElementwisePhysicalAdvPseudoPhysicalFootFinderMem**_](classpolar__foot__finder__details_1_1ElementwisePhysicalAdvPseudoPhysicalFootFinderMem.md) _._
```C++
inline polar_foot_finder_details::ElementwisePhysicalAdvPseudoPhysicalFootFinderMem::ElementwisePhysicalAdvPseudoPhysicalFootFinderMem (
    RThetaAdvectionEvaluator const & evaluator_advection_field,
    LogicalToPhysicalMapping const & logical_to_physical,
    TimeStepperBuilder const & time_stepper_builder,
    AdvecCoefFieldMem && advection_field_coefs,
    Coord< X_pC , Y_pC > coord_centre,
    IdxRange< GridTheta > idx_range_theta,
    double epsilon=1e-12
) 
```





**Parameters:**


* `evaluator_advection_field` The evaluator for the spline representation of the advection field. 
* `logical_to_physical` The mapping from the logical domain to the physical domain, used to construct the pseudo-physical-to-advection mapping. 
* `time_stepper_builder` The factory used to preallocate the time integration method. 
* `advection_field_coefs` The spline coefficients of the advection field. Ownership is transferred in. 
* `coord_centre` The coordinate of the polar centre in the pseudo-physical domain, used to handle the degenerate point at \(r = 0\). 
* `idx_range_theta` The poloidal index range, used to wrap the angular coordinate into the periodic domain after each time step. 
* `epsilon` \(\varepsilon\) parameter used for the linearisation of the advection field around the central point. 




        

<hr>



### function operator() 

_Create an_ [_**ElementwisePhysicalAdvPseudoPhysicalFootFinder**_](classpolar__foot__finder__details_1_1ElementwisePhysicalAdvPseudoPhysicalFootFinder.md) _for the given time step._
```C++
inline GPUCompat polar_foot_finder_details::ElementwisePhysicalAdvPseudoPhysicalFootFinderMem::operator() (
    double dt
) 
```



Returns a non-owning [**ElementwisePhysicalAdvPseudoPhysicalFootFinder**](classpolar__foot__finder__details_1_1ElementwisePhysicalAdvPseudoPhysicalFootFinder.md) that holds views of the stored spline coefficients and is configured for time step \(dt\). The returned object can be copied to and called from the device.




**Parameters:**


* `dt` The time step for the characteristic equation.



**Returns:**

A view-based [**ElementwisePhysicalAdvPseudoPhysicalFootFinder**](classpolar__foot__finder__details_1_1ElementwisePhysicalAdvPseudoPhysicalFootFinder.md) for the given time step. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/advection/polar_foot_finders/physical_advection_pseudo_physical_foot_finder.hpp`

