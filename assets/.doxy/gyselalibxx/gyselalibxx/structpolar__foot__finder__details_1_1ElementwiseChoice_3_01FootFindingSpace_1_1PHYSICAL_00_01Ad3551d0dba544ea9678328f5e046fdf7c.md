

# Struct polar\_foot\_finder\_details::ElementwiseChoice&lt; FootFindingSpace::PHYSICAL, AdvectionFieldSpace::LOGICAL, GridR, GridTheta, IdxRangeOperator, RThetaAdvectionEvaluator, AdvecCoefField, TimeStepperBuilder, LogicalToPhysicalMapping &gt;

**template &lt;class [**GridR**](structGridR.md), class [**GridTheta**](structGridTheta.md), class IdxRangeOperator, class RThetaAdvectionEvaluator, class AdvecCoefField, class TimeStepperBuilder, concepts::Mapping LogicalToPhysicalMapping&gt;**



[**ClassList**](annotated.md) **>** [**polar\_foot\_finder\_details**](namespacepolar__foot__finder__details.md) **>** [**ElementwiseChoice&lt; FootFindingSpace::PHYSICAL, AdvectionFieldSpace::LOGICAL, GridR, GridTheta, IdxRangeOperator, RThetaAdvectionEvaluator, AdvecCoefField, TimeStepperBuilder, LogicalToPhysicalMapping &gt;**](structpolar__foot__finder__details_1_1ElementwiseChoice_3_01FootFindingSpace_1_1PHYSICAL_00_01Ad3551d0dba544ea9678328f5e046fdf7c.md)



_Selects_ [_**ElementwiseLogicalAdvPseudoPhysFootFinderMem**_](classpolar__foot__finder__details_1_1ElementwiseLogicalAdvPseudoPhysFootFinderMem.md) _(using the physical mapping directly) for logical advection with foot-finding in physical space._

* `#include <logical_advection_pseudo_physical_foot_finder.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef [**ElementwiseLogicalAdvPseudoPhysFootFinderMem**](classpolar__foot__finder__details_1_1ElementwiseLogicalAdvPseudoPhysFootFinderMem.md)&lt; [**GridR**](structGridR.md), [**GridTheta**](structGridTheta.md), IdxRangeOperator, RThetaAdvectionEvaluator, AdvecCoefField, TimeStepperBuilder, LogicalToPhysicalMapping &gt; | [**type**](#typedef-type)  <br>_The selected elementwise operator type._  |
















































## Public Types Documentation




### typedef type 

_The selected elementwise operator type._ 
```C++
using polar_foot_finder_details::ElementwiseChoice< FootFindingSpace::PHYSICAL, AdvectionFieldSpace::LOGICAL, GridR, GridTheta, IdxRangeOperator, RThetaAdvectionEvaluator, AdvecCoefField, TimeStepperBuilder, LogicalToPhysicalMapping >::type =  ElementwiseLogicalAdvPseudoPhysFootFinderMem< GridR, GridTheta, IdxRangeOperator, RThetaAdvectionEvaluator, AdvecCoefField, TimeStepperBuilder, LogicalToPhysicalMapping>;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/advection/polar_foot_finders/logical_advection_pseudo_physical_foot_finder.hpp`

