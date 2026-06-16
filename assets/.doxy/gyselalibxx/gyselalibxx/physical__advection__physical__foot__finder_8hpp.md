

# File physical\_advection\_physical\_foot\_finder.hpp



[**FileList**](files.md) **>** [**advection**](dir_b90fde0f10c67a9aef841a6e6700f1f6.md) **>** [**polar\_foot\_finders**](dir_71ecba75f1675b0dc7d7b6463874ff90.md) **>** [**physical\_advection\_physical\_foot\_finder.hpp**](physical__advection__physical__foot__finder_8hpp.md)

[Go to the source code of this file](physical__advection__physical__foot__finder_8hpp_source.md)



* `#include <ddc/ddc.hpp>`
* `#include "ddc_aliases.hpp"`
* `#include "ddc_helper.hpp"`
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
| struct | [**ElementwiseChoice&lt; FootFindingSpace::PHYSICAL, AdvectionFieldSpace::PHYSICAL, GridR, GridTheta, IdxRangeOperator, RThetaAdvectionEvaluator, AdvecCoefField, TimeStepperBuilder, LogicalToPhysicalMapping &gt;**](structpolar__foot__finder__details_1_1ElementwiseChoice_3_01FootFindingSpace_1_1PHYSICAL_00_01Ade857839a3fb92baaf6bd919a12083d58.md) &lt;class [**GridR**](structGridR.md), class [**GridTheta**](structGridTheta.md), class IdxRangeOperator, class RThetaAdvectionEvaluator, class AdvecCoefField, class TimeStepperBuilder, LogicalToPhysicalMapping&gt;<br>_Selects_ [_**ElementwisePhysicalAdvPhysicalFootFinderMem**_](classpolar__foot__finder__details_1_1ElementwisePhysicalAdvPhysicalFootFinderMem.md) _for physical advection with foot-finding in physical space._ |
| class | [**ElementwisePhysicalAdvPhysicalFootFinder**](classpolar__foot__finder__details_1_1ElementwisePhysicalAdvPhysicalFootFinder.md) &lt;class [**GridR**](structGridR.md), class [**GridTheta**](structGridTheta.md), class IdxRangeOperator, class RThetaAdvectionEvaluator, class AdvecCoefField, class TimeStepper, class LogicalToPhysicalMapping&gt;<br>_GPU-callable functor that finds characteristic feet for physical-space advection with foot-finding in physical_ \((x, y)\) _space._ |
| class | [**ElementwisePhysicalAdvPhysicalFootFinderMem**](classpolar__foot__finder__details_1_1ElementwisePhysicalAdvPhysicalFootFinderMem.md) &lt;class [**GridR**](structGridR.md), class [**GridTheta**](structGridTheta.md), class IdxRangeOperator, class RThetaAdvectionEvaluator, class AdvecCoefFieldMem, class TimeStepperBuilder, class LogicalToPhysicalMapping&gt;<br>_Owning manager for_ [_**ElementwisePhysicalAdvPhysicalFootFinder**_](classpolar__foot__finder__details_1_1ElementwisePhysicalAdvPhysicalFootFinder.md) _._ |



















































------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/advection/polar_foot_finders/physical_advection_physical_foot_finder.hpp`

