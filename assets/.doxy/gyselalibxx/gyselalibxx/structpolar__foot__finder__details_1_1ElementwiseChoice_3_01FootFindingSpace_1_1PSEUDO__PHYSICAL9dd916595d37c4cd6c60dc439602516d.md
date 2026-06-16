

# Struct polar\_foot\_finder\_details::ElementwiseChoice&lt; FootFindingSpace::PSEUDO\_PHYSICAL, AdvectionFieldSpace::LOGICAL, GridR, GridTheta, IdxRangeOperator, RThetaAdvectionEvaluator, AdvecCoefField, TimeStepperBuilder, LogicalToPhysicalMapping &gt;

**template &lt;class [**GridR**](structGridR.md), class [**GridTheta**](structGridTheta.md), class IdxRangeOperator, class RThetaAdvectionEvaluator, class AdvecCoefField, class TimeStepperBuilder, concepts::Mapping LogicalToPhysicalMapping&gt;**



[**ClassList**](annotated.md) **>** [**polar\_foot\_finder\_details**](namespacepolar__foot__finder__details.md) **>** [**ElementwiseChoice&lt; FootFindingSpace::PSEUDO\_PHYSICAL, AdvectionFieldSpace::LOGICAL, GridR, GridTheta, IdxRangeOperator, RThetaAdvectionEvaluator, AdvecCoefField, TimeStepperBuilder, LogicalToPhysicalMapping &gt;**](structpolar__foot__finder__details_1_1ElementwiseChoice_3_01FootFindingSpace_1_1PSEUDO__PHYSICAL9dd916595d37c4cd6c60dc439602516d.md)



_Selects_ [_**ElementwiseLogicalAdvPseudoPhysFootFinderMem**_](classpolar__foot__finder__details_1_1ElementwiseLogicalAdvPseudoPhysFootFinderMem.md) _(circular pseudo-physical mapping) for logical advection with foot-finding in pseudo-physical space._

* `#include <logical_advection_pseudo_physical_foot_finder.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef [**ElementwiseLogicalAdvPseudoPhysFootFinderMem**](classpolar__foot__finder__details_1_1ElementwiseLogicalAdvPseudoPhysFootFinderMem.md)&lt; [**GridR**](structGridR.md), [**GridTheta**](structGridTheta.md), IdxRangeOperator, RThetaAdvectionEvaluator, AdvecCoefField, TimeStepperBuilder, [**CircularToCartesian**](classCircularToCartesian.md)&lt; typename GridR::continuous\_dimension\_type, typename GridTheta::continuous\_dimension\_type, [**X\_pC**](structX__pC.md), [**Y\_pC**](structY__pC.md) &gt; &gt; | [**type**](#typedef-type)  <br>_The selected elementwise operator type._  |
















































## Public Types Documentation




### typedef type 

_The selected elementwise operator type._ 
```C++
using polar_foot_finder_details::ElementwiseChoice< FootFindingSpace::PSEUDO_PHYSICAL, AdvectionFieldSpace::LOGICAL, GridR, GridTheta, IdxRangeOperator, RThetaAdvectionEvaluator, AdvecCoefField, TimeStepperBuilder, LogicalToPhysicalMapping >::type =  ElementwiseLogicalAdvPseudoPhysFootFinderMem< GridR, GridTheta, IdxRangeOperator, RThetaAdvectionEvaluator, AdvecCoefField, TimeStepperBuilder, CircularToCartesian< typename GridR::continuous_dimension_type, typename GridTheta::continuous_dimension_type, X_pC, Y_pC> >;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/advection/polar_foot_finders/logical_advection_pseudo_physical_foot_finder.hpp`

