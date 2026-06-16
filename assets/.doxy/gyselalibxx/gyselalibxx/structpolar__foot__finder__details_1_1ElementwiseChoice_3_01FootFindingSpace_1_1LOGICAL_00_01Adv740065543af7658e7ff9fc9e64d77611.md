

# Struct polar\_foot\_finder\_details::ElementwiseChoice&lt; FootFindingSpace::LOGICAL, AdvectionFieldSpace::LOGICAL, GridR, GridTheta, IdxRangeOperator, RThetaAdvectionEvaluator, AdvecCoefField, TimeStepperBuilder, LogicalToPhysicalMapping &gt;

**template &lt;class [**GridR**](structGridR.md), class [**GridTheta**](structGridTheta.md), class IdxRangeOperator, class RThetaAdvectionEvaluator, class AdvecCoefField, class TimeStepperBuilder, concepts::Mapping LogicalToPhysicalMapping&gt;**



[**ClassList**](annotated.md) **>** [**polar\_foot\_finder\_details**](namespacepolar__foot__finder__details.md) **>** [**ElementwiseChoice&lt; FootFindingSpace::LOGICAL, AdvectionFieldSpace::LOGICAL, GridR, GridTheta, IdxRangeOperator, RThetaAdvectionEvaluator, AdvecCoefField, TimeStepperBuilder, LogicalToPhysicalMapping &gt;**](structpolar__foot__finder__details_1_1ElementwiseChoice_3_01FootFindingSpace_1_1LOGICAL_00_01Adv740065543af7658e7ff9fc9e64d77611.md)



_Selects_ [_**ElementwiseLogicalAdvLogicalFootFinderMem**_](classpolar__foot__finder__details_1_1ElementwiseLogicalAdvLogicalFootFinderMem.md) _for logical advection with foot-finding in logical space._

* `#include <logical_advection_logical_foot_finder.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef [**ElementwiseLogicalAdvLogicalFootFinderMem**](classpolar__foot__finder__details_1_1ElementwiseLogicalAdvLogicalFootFinderMem.md)&lt; [**GridR**](structGridR.md), [**GridTheta**](structGridTheta.md), IdxRangeOperator, RThetaAdvectionEvaluator, AdvecCoefField, TimeStepperBuilder &gt; | [**type**](#typedef-type)  <br>_The selected elementwise operator type._  |
















































## Public Types Documentation




### typedef type 

_The selected elementwise operator type._ 
```C++
using polar_foot_finder_details::ElementwiseChoice< FootFindingSpace::LOGICAL, AdvectionFieldSpace::LOGICAL, GridR, GridTheta, IdxRangeOperator, RThetaAdvectionEvaluator, AdvecCoefField, TimeStepperBuilder, LogicalToPhysicalMapping >::type =  ElementwiseLogicalAdvLogicalFootFinderMem< GridR, GridTheta, IdxRangeOperator, RThetaAdvectionEvaluator, AdvecCoefField, TimeStepperBuilder>;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/advection/polar_foot_finders/logical_advection_logical_foot_finder.hpp`

