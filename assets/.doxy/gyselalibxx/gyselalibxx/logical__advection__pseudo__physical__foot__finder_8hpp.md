

# File logical\_advection\_pseudo\_physical\_foot\_finder.hpp



[**FileList**](files.md) **>** [**advection**](dir_b90fde0f10c67a9aef841a6e6700f1f6.md) **>** [**polar\_foot\_finders**](dir_71ecba75f1675b0dc7d7b6463874ff90.md) **>** [**logical\_advection\_pseudo\_physical\_foot\_finder.hpp**](logical__advection__pseudo__physical__foot__finder_8hpp.md)

[Go to the source code of this file](logical__advection__pseudo__physical__foot__finder_8hpp_source.md)



* `#include <ddc/ddc.hpp>`
* `#include "circular_to_cartesian.hpp"`
* `#include "ddc_aliases.hpp"`
* `#include "ddc_helper.hpp"`
* `#include "geometry_pseudo_cartesian.hpp"`
* `#include "l_norm_tools.hpp"`
* `#include "tensor.hpp"`
* `#include "type_seq_tools.hpp"`
* `#include "vector_field.hpp"`













## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**polar\_foot\_finder\_details**](namespacepolar__foot__finder__details.md) <br> |


## Classes

| Type | Name |
| ---: | :--- |
| struct | [**ElementwiseChoice&lt; FootFindingSpace::PHYSICAL, AdvectionFieldSpace::LOGICAL, GridR, GridTheta, IdxRangeOperator, RThetaAdvectionEvaluator, AdvecCoefField, TimeStepperBuilder, LogicalToPhysicalMapping &gt;**](structpolar__foot__finder__details_1_1ElementwiseChoice_3_01FootFindingSpace_1_1PHYSICAL_00_01Ad3551d0dba544ea9678328f5e046fdf7c.md) &lt;class [**GridR**](structGridR.md), class [**GridTheta**](structGridTheta.md), class IdxRangeOperator, class RThetaAdvectionEvaluator, class AdvecCoefField, class TimeStepperBuilder, LogicalToPhysicalMapping&gt;<br>_Selects_ [_**ElementwiseLogicalAdvPseudoPhysFootFinderMem**_](classpolar__foot__finder__details_1_1ElementwiseLogicalAdvPseudoPhysFootFinderMem.md) _(using the physical mapping directly) for logical advection with foot-finding in physical space._ |
| struct | [**ElementwiseChoice&lt; FootFindingSpace::PSEUDO\_PHYSICAL, AdvectionFieldSpace::LOGICAL, GridR, GridTheta, IdxRangeOperator, RThetaAdvectionEvaluator, AdvecCoefField, TimeStepperBuilder, LogicalToPhysicalMapping &gt;**](structpolar__foot__finder__details_1_1ElementwiseChoice_3_01FootFindingSpace_1_1PSEUDO__PHYSICAL9dd916595d37c4cd6c60dc439602516d.md) &lt;class [**GridR**](structGridR.md), class [**GridTheta**](structGridTheta.md), class IdxRangeOperator, class RThetaAdvectionEvaluator, class AdvecCoefField, class TimeStepperBuilder, LogicalToPhysicalMapping&gt;<br>_Selects_ [_**ElementwiseLogicalAdvPseudoPhysFootFinderMem**_](classpolar__foot__finder__details_1_1ElementwiseLogicalAdvPseudoPhysFootFinderMem.md) _(circular pseudo-physical mapping) for logical advection with foot-finding in pseudo-physical space._ |
| class | [**ElementwiseLogicalAdvPseudoPhysFootFinder**](classpolar__foot__finder__details_1_1ElementwiseLogicalAdvPseudoPhysFootFinder.md) &lt;class [**GridR**](structGridR.md), class [**GridTheta**](structGridTheta.md), class IdxRangeOperator, class RThetaAdvectionEvaluator, class AdvecCoefField, class TimeStepper, class LogicalToPseudoPhysicalMapping&gt;<br>_GPU-callable functor that finds characteristic feet for logical-space advection with foot-finding in pseudo-physical space._  |
| class | [**ElementwiseLogicalAdvPseudoPhysFootFinderMem**](classpolar__foot__finder__details_1_1ElementwiseLogicalAdvPseudoPhysFootFinderMem.md) &lt;class [**GridR**](structGridR.md), class [**GridTheta**](structGridTheta.md), class IdxRangeOperator, class RThetaAdvectionEvaluator, class AdvecCoefFieldMem, class TimeStepperBuilder, class LogicalToPseudoPhysicalMapping&gt;<br>_Owning manager for_ [_**ElementwiseLogicalAdvPseudoPhysFootFinder**_](classpolar__foot__finder__details_1_1ElementwiseLogicalAdvPseudoPhysFootFinder.md) _._ |



















































------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/advection/polar_foot_finders/logical_advection_pseudo_physical_foot_finder.hpp`

