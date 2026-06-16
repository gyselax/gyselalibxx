

# Namespace polar\_foot\_finder\_details



[**Namespace List**](namespaces.md) **>** [**polar\_foot\_finder\_details**](namespacepolar__foot__finder__details.md)




















## Classes

| Type | Name |
| ---: | :--- |
| class | [**ElementwiseChoice**](classpolar__foot__finder__details_1_1ElementwiseChoice.md) &lt;FFSpace, AFSpace, class [**GridR**](structGridR.md), class [**GridTheta**](structGridTheta.md), class IdxRangeOperator, class RThetaAdvectionEvaluator, class AdvecCoefField, class TimeStepperBuilder, LogicalToPhysicalMapping&gt;<br> |
| struct | [**ElementwiseChoice&lt; FootFindingSpace::LOGICAL, AdvectionFieldSpace::LOGICAL, GridR, GridTheta, IdxRangeOperator, RThetaAdvectionEvaluator, AdvecCoefField, TimeStepperBuilder, LogicalToPhysicalMapping &gt;**](structpolar__foot__finder__details_1_1ElementwiseChoice_3_01FootFindingSpace_1_1LOGICAL_00_01Adv740065543af7658e7ff9fc9e64d77611.md) &lt;class [**GridR**](structGridR.md), class [**GridTheta**](structGridTheta.md), class IdxRangeOperator, class RThetaAdvectionEvaluator, class AdvecCoefField, class TimeStepperBuilder, LogicalToPhysicalMapping&gt;<br>_Selects_ [_**ElementwiseLogicalAdvLogicalFootFinderMem**_](classpolar__foot__finder__details_1_1ElementwiseLogicalAdvLogicalFootFinderMem.md) _for logical advection with foot-finding in logical space._ |
| struct | [**ElementwiseChoice&lt; FootFindingSpace::PHYSICAL, AdvectionFieldSpace::LOGICAL, GridR, GridTheta, IdxRangeOperator, RThetaAdvectionEvaluator, AdvecCoefField, TimeStepperBuilder, LogicalToPhysicalMapping &gt;**](structpolar__foot__finder__details_1_1ElementwiseChoice_3_01FootFindingSpace_1_1PHYSICAL_00_01Ad3551d0dba544ea9678328f5e046fdf7c.md) &lt;class [**GridR**](structGridR.md), class [**GridTheta**](structGridTheta.md), class IdxRangeOperator, class RThetaAdvectionEvaluator, class AdvecCoefField, class TimeStepperBuilder, LogicalToPhysicalMapping&gt;<br>_Selects_ [_**ElementwiseLogicalAdvPseudoPhysFootFinderMem**_](classpolar__foot__finder__details_1_1ElementwiseLogicalAdvPseudoPhysFootFinderMem.md) _(using the physical mapping directly) for logical advection with foot-finding in physical space._ |
| struct | [**ElementwiseChoice&lt; FootFindingSpace::PHYSICAL, AdvectionFieldSpace::PHYSICAL, GridR, GridTheta, IdxRangeOperator, RThetaAdvectionEvaluator, AdvecCoefField, TimeStepperBuilder, LogicalToPhysicalMapping &gt;**](structpolar__foot__finder__details_1_1ElementwiseChoice_3_01FootFindingSpace_1_1PHYSICAL_00_01Ade857839a3fb92baaf6bd919a12083d58.md) &lt;class [**GridR**](structGridR.md), class [**GridTheta**](structGridTheta.md), class IdxRangeOperator, class RThetaAdvectionEvaluator, class AdvecCoefField, class TimeStepperBuilder, LogicalToPhysicalMapping&gt;<br>_Selects_ [_**ElementwisePhysicalAdvPhysicalFootFinderMem**_](classpolar__foot__finder__details_1_1ElementwisePhysicalAdvPhysicalFootFinderMem.md) _for physical advection with foot-finding in physical space._ |
| struct | [**ElementwiseChoice&lt; FootFindingSpace::PSEUDO\_PHYSICAL, AdvectionFieldSpace::LOGICAL, GridR, GridTheta, IdxRangeOperator, RThetaAdvectionEvaluator, AdvecCoefField, TimeStepperBuilder, LogicalToPhysicalMapping &gt;**](structpolar__foot__finder__details_1_1ElementwiseChoice_3_01FootFindingSpace_1_1PSEUDO__PHYSICAL9dd916595d37c4cd6c60dc439602516d.md) &lt;class [**GridR**](structGridR.md), class [**GridTheta**](structGridTheta.md), class IdxRangeOperator, class RThetaAdvectionEvaluator, class AdvecCoefField, class TimeStepperBuilder, LogicalToPhysicalMapping&gt;<br>_Selects_ [_**ElementwiseLogicalAdvPseudoPhysFootFinderMem**_](classpolar__foot__finder__details_1_1ElementwiseLogicalAdvPseudoPhysFootFinderMem.md) _(circular pseudo-physical mapping) for logical advection with foot-finding in pseudo-physical space._ |
| struct | [**ElementwiseChoice&lt; FootFindingSpace::PSEUDO\_PHYSICAL, AdvectionFieldSpace::PHYSICAL, GridR, GridTheta, IdxRangeOperator, RThetaAdvectionEvaluator, AdvecCoefField, TimeStepperBuilder, LogicalToPhysicalMapping &gt;**](structpolar__foot__finder__details_1_1ElementwiseChoice_3_01FootFindingSpace_1_1PSEUDO__PHYSICAL0e03af03e715868452e393a1a59b4674.md) &lt;class [**GridR**](structGridR.md), class [**GridTheta**](structGridTheta.md), class IdxRangeOperator, class RThetaAdvectionEvaluator, class AdvecCoefField, class TimeStepperBuilder, LogicalToPhysicalMapping&gt;<br>_Selects_ [_**ElementwisePhysicalAdvPseudoPhysicalFootFinderMem**_](classpolar__foot__finder__details_1_1ElementwisePhysicalAdvPseudoPhysicalFootFinderMem.md) _for physical advection with foot-finding in pseudo-physical space._ |
| class | [**ElementwiseLogicalAdvLogicalFootFinder**](classpolar__foot__finder__details_1_1ElementwiseLogicalAdvLogicalFootFinder.md) &lt;class [**GridR**](structGridR.md), class [**GridTheta**](structGridTheta.md), class IdxRangeOperator, class RThetaAdvectionEvaluator, class AdvecCoefField, class TimeStepper&gt;<br>_GPU-callable functor that finds characteristic feet for logical-space advection with foot-finding in logical_ \((r, \theta)\) _space._ |
| class | [**ElementwiseLogicalAdvLogicalFootFinderMem**](classpolar__foot__finder__details_1_1ElementwiseLogicalAdvLogicalFootFinderMem.md) &lt;class [**GridR**](structGridR.md), class [**GridTheta**](structGridTheta.md), class IdxRangeOperator, class RThetaAdvectionEvaluator, class AdvecCoefFieldMem, class TimeStepperBuilder&gt;<br>_Owning manager for_ [_**ElementwiseLogicalAdvLogicalFootFinder**_](classpolar__foot__finder__details_1_1ElementwiseLogicalAdvLogicalFootFinder.md) _._ |
| class | [**ElementwiseLogicalAdvPseudoPhysFootFinder**](classpolar__foot__finder__details_1_1ElementwiseLogicalAdvPseudoPhysFootFinder.md) &lt;class [**GridR**](structGridR.md), class [**GridTheta**](structGridTheta.md), class IdxRangeOperator, class RThetaAdvectionEvaluator, class AdvecCoefField, class TimeStepper, class LogicalToPseudoPhysicalMapping&gt;<br>_GPU-callable functor that finds characteristic feet for logical-space advection with foot-finding in pseudo-physical space._  |
| class | [**ElementwiseLogicalAdvPseudoPhysFootFinderMem**](classpolar__foot__finder__details_1_1ElementwiseLogicalAdvPseudoPhysFootFinderMem.md) &lt;class [**GridR**](structGridR.md), class [**GridTheta**](structGridTheta.md), class IdxRangeOperator, class RThetaAdvectionEvaluator, class AdvecCoefFieldMem, class TimeStepperBuilder, class LogicalToPseudoPhysicalMapping&gt;<br>_Owning manager for_ [_**ElementwiseLogicalAdvPseudoPhysFootFinder**_](classpolar__foot__finder__details_1_1ElementwiseLogicalAdvPseudoPhysFootFinder.md) _._ |
| class | [**ElementwisePhysicalAdvPhysicalFootFinder**](classpolar__foot__finder__details_1_1ElementwisePhysicalAdvPhysicalFootFinder.md) &lt;class [**GridR**](structGridR.md), class [**GridTheta**](structGridTheta.md), class IdxRangeOperator, class RThetaAdvectionEvaluator, class AdvecCoefField, class TimeStepper, class LogicalToPhysicalMapping&gt;<br>_GPU-callable functor that finds characteristic feet for physical-space advection with foot-finding in physical_ \((x, y)\) _space._ |
| class | [**ElementwisePhysicalAdvPhysicalFootFinderMem**](classpolar__foot__finder__details_1_1ElementwisePhysicalAdvPhysicalFootFinderMem.md) &lt;class [**GridR**](structGridR.md), class [**GridTheta**](structGridTheta.md), class IdxRangeOperator, class RThetaAdvectionEvaluator, class AdvecCoefFieldMem, class TimeStepperBuilder, class LogicalToPhysicalMapping&gt;<br>_Owning manager for_ [_**ElementwisePhysicalAdvPhysicalFootFinder**_](classpolar__foot__finder__details_1_1ElementwisePhysicalAdvPhysicalFootFinder.md) _._ |
| class | [**ElementwisePhysicalAdvPseudoPhysicalFootFinder**](classpolar__foot__finder__details_1_1ElementwisePhysicalAdvPseudoPhysicalFootFinder.md) &lt;class [**GridR**](structGridR.md), class [**GridTheta**](structGridTheta.md), class IdxRangeOperator, class RThetaAdvectionEvaluator, class PseudoPhysicalToAdvectionMapping, class PseudoPhysicalToLogicalMapping, class LogicalToPseudoPhysicalMapping, class AdvecCoefField, class TimeStepper&gt;<br>_GPU-callable functor that finds characteristic feet for physical-space advection with foot-finding in pseudo-physical_ \((X_{pC}, Y_{pC})\) _space._ |
| class | [**ElementwisePhysicalAdvPseudoPhysicalFootFinderMem**](classpolar__foot__finder__details_1_1ElementwisePhysicalAdvPseudoPhysicalFootFinderMem.md) &lt;class [**GridR**](structGridR.md), class [**GridTheta**](structGridTheta.md), class IdxRangeOperator, class RThetaAdvectionEvaluator, class AdvecCoefFieldMem, class TimeStepperBuilder, LogicalToPhysicalMapping&gt;<br>_Owning manager for_ [_**ElementwisePhysicalAdvPseudoPhysicalFootFinder**_](classpolar__foot__finder__details_1_1ElementwisePhysicalAdvPseudoPhysicalFootFinder.md) _._ |



















































------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/advection/polar_foot_finders/elementwise_choice.hpp`

