

# Struct polar\_foot\_finder\_details::ElementwiseChoice&lt; FootFindingSpace::PHYSICAL, AdvectionFieldSpace::PHYSICAL, GridR, GridTheta, IdxRangeOperator, RThetaAdvectionEvaluator, AdvecCoefField, TimeStepperBuilder, LogicalToPhysicalMapping &gt;

**template &lt;class [**GridR**](structGridR.md), class [**GridTheta**](structGridTheta.md), class IdxRangeOperator, class RThetaAdvectionEvaluator, class AdvecCoefField, class TimeStepperBuilder, concepts::Mapping LogicalToPhysicalMapping&gt;**



[**ClassList**](annotated.md) **>** [**polar\_foot\_finder\_details**](namespacepolar__foot__finder__details.md) **>** [**ElementwiseChoice&lt; FootFindingSpace::PHYSICAL, AdvectionFieldSpace::PHYSICAL, GridR, GridTheta, IdxRangeOperator, RThetaAdvectionEvaluator, AdvecCoefField, TimeStepperBuilder, LogicalToPhysicalMapping &gt;**](structpolar__foot__finder__details_1_1ElementwiseChoice_3_01FootFindingSpace_1_1PHYSICAL_00_01Ade857839a3fb92baaf6bd919a12083d58.md)



_Selects_ [_**ElementwisePhysicalAdvPhysicalFootFinderMem**_](classpolar__foot__finder__details_1_1ElementwisePhysicalAdvPhysicalFootFinderMem.md) _for physical advection with foot-finding in physical space._

* `#include <physical_advection_physical_foot_finder.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef [**ElementwisePhysicalAdvPhysicalFootFinderMem**](classpolar__foot__finder__details_1_1ElementwisePhysicalAdvPhysicalFootFinderMem.md)&lt; [**GridR**](structGridR.md), [**GridTheta**](structGridTheta.md), IdxRangeOperator, RThetaAdvectionEvaluator, AdvecCoefField, TimeStepperBuilder, LogicalToPhysicalMapping &gt; | [**type**](#typedef-type)  <br>_The selected elementwise operator type._  |
















































## Public Types Documentation




### typedef type 

_The selected elementwise operator type._ 
```C++
using polar_foot_finder_details::ElementwiseChoice< FootFindingSpace::PHYSICAL, AdvectionFieldSpace::PHYSICAL, GridR, GridTheta, IdxRangeOperator, RThetaAdvectionEvaluator, AdvecCoefField, TimeStepperBuilder, LogicalToPhysicalMapping >::type =  ElementwisePhysicalAdvPhysicalFootFinderMem< GridR, GridTheta, IdxRangeOperator, RThetaAdvectionEvaluator, AdvecCoefField, TimeStepperBuilder, LogicalToPhysicalMapping>;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/advection/polar_foot_finders/physical_advection_physical_foot_finder.hpp`

