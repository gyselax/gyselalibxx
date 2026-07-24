

# File logical\_advection\_logical\_foot\_finder.hpp



[**FileList**](files.md) **>** [**advection**](dir_b90fde0f10c67a9aef841a6e6700f1f6.md) **>** [**polar\_foot\_finders**](dir_71ecba75f1675b0dc7d7b6463874ff90.md) **>** [**logical\_advection\_logical\_foot\_finder.hpp**](logical__advection__logical__foot__finder_8hpp.md)

[Go to the source code of this file](logical__advection__logical__foot__finder_8hpp_source.md)



* `#include <ddc/ddc.hpp>`
* `#include "coord_transformation_tools.hpp"`
* `#include "ddc_aliases.hpp"`
* `#include "ddc_helper.hpp"`
* `#include "elementwise_choice.hpp"`
* `#include "l_norm_tools.hpp"`
* `#include "tensor.hpp"`
* `#include "vector_field.hpp"`













## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**polar\_foot\_finder\_details**](namespacepolar__foot__finder__details.md) <br> |


## Classes

| Type | Name |
| ---: | :--- |
| struct | [**ElementwiseChoice&lt; FootFindingSpace::LOGICAL, AdvectionFieldSpace::LOGICAL, GridR, GridTheta, IdxRangeOperator, RThetaAdvectionEvaluator, AdvecCoefField, TimeStepperBuilder, LogicalToPhysicalMapping &gt;**](structpolar__foot__finder__details_1_1ElementwiseChoice_3_01FootFindingSpace_1_1LOGICAL_00_01Adv740065543af7658e7ff9fc9e64d77611.md) &lt;class [**GridR**](structGridR.md), class [**GridTheta**](structGridTheta.md), class IdxRangeOperator, class RThetaAdvectionEvaluator, class AdvecCoefField, class TimeStepperBuilder, LogicalToPhysicalMapping&gt;<br>_Selects_ [_**ElementwiseLogicalAdvLogicalFootFinderMem**_](classpolar__foot__finder__details_1_1ElementwiseLogicalAdvLogicalFootFinderMem.md) _for logical advection with foot-finding in logical space._ |
| class | [**ElementwiseLogicalAdvLogicalFootFinder**](classpolar__foot__finder__details_1_1ElementwiseLogicalAdvLogicalFootFinder.md) &lt;class [**GridR**](structGridR.md), class [**GridTheta**](structGridTheta.md), class IdxRangeOperator, class RThetaAdvectionEvaluator, class AdvecCoefField, class TimeStepper&gt;<br>_GPU-callable functor that finds characteristic feet for logical-space advection with foot-finding in logical_ \((r, \theta)\) _space._ |
| class | [**ElementwiseLogicalAdvLogicalFootFinderMem**](classpolar__foot__finder__details_1_1ElementwiseLogicalAdvLogicalFootFinderMem.md) &lt;class [**GridR**](structGridR.md), class [**GridTheta**](structGridTheta.md), class IdxRangeOperator, class RThetaAdvectionEvaluator, class AdvecCoefFieldMem, class TimeStepperBuilder&gt;<br>_Owning manager for_ [_**ElementwiseLogicalAdvLogicalFootFinder**_](classpolar__foot__finder__details_1_1ElementwiseLogicalAdvLogicalFootFinder.md) _._ |



















































------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/advection/polar_foot_finders/logical_advection_logical_foot_finder.hpp`

