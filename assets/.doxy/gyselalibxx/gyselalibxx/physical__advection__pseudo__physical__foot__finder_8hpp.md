

# File physical\_advection\_pseudo\_physical\_foot\_finder.hpp



[**FileList**](files.md) **>** [**advection**](dir_b90fde0f10c67a9aef841a6e6700f1f6.md) **>** [**polar\_foot\_finders**](dir_71ecba75f1675b0dc7d7b6463874ff90.md) **>** [**physical\_advection\_pseudo\_physical\_foot\_finder.hpp**](physical__advection__pseudo__physical__foot__finder_8hpp.md)

[Go to the source code of this file](physical__advection__pseudo__physical__foot__finder_8hpp_source.md)



* `#include <ddc/ddc.hpp>`
* `#include "combined_mapping.hpp"`
* `#include "ddc_aliases.hpp"`
* `#include "ddc_helper.hpp"`
* `#include "geometry_pseudo_cartesian.hpp"`
* `#include "l_norm_tools.hpp"`
* `#include "tensor.hpp"`
* `#include "type_seq_tools.hpp"`
* `#include "vector_field.hpp"`
* `#include "vector_mapper.hpp"`













## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**polar\_foot\_finder\_details**](namespacepolar__foot__finder__details.md) <br> |


## Classes

| Type | Name |
| ---: | :--- |
| struct | [**ElementwiseChoice&lt; FootFindingSpace::PSEUDO\_PHYSICAL, AdvectionFieldSpace::PHYSICAL, GridR, GridTheta, IdxRangeOperator, RThetaAdvectionEvaluator, AdvecCoefField, TimeStepperBuilder, LogicalToPhysicalMapping &gt;**](structpolar__foot__finder__details_1_1ElementwiseChoice_3_01FootFindingSpace_1_1PSEUDO__PHYSICAL0e03af03e715868452e393a1a59b4674.md) &lt;class [**GridR**](structGridR.md), class [**GridTheta**](structGridTheta.md), class IdxRangeOperator, class RThetaAdvectionEvaluator, class AdvecCoefField, class TimeStepperBuilder, LogicalToPhysicalMapping&gt;<br>_Selects_ [_**ElementwisePhysicalAdvPseudoPhysicalFootFinderMem**_](classpolar__foot__finder__details_1_1ElementwisePhysicalAdvPseudoPhysicalFootFinderMem.md) _for physical advection with foot-finding in pseudo-physical space._ |
| class | [**ElementwisePhysicalAdvPseudoPhysicalFootFinder**](classpolar__foot__finder__details_1_1ElementwisePhysicalAdvPseudoPhysicalFootFinder.md) &lt;class [**GridR**](structGridR.md), class [**GridTheta**](structGridTheta.md), class IdxRangeOperator, class RThetaAdvectionEvaluator, class PseudoPhysicalToAdvectionMapping, class PseudoPhysicalToLogicalMapping, class LogicalToPseudoPhysicalMapping, class AdvecCoefField, class TimeStepper&gt;<br>_GPU-callable functor that finds characteristic feet for physical-space advection with foot-finding in pseudo-physical_ \((X_{pC}, Y_{pC})\) _space._ |
| class | [**ElementwisePhysicalAdvPseudoPhysicalFootFinderMem**](classpolar__foot__finder__details_1_1ElementwisePhysicalAdvPseudoPhysicalFootFinderMem.md) &lt;class [**GridR**](structGridR.md), class [**GridTheta**](structGridTheta.md), class IdxRangeOperator, class RThetaAdvectionEvaluator, class AdvecCoefFieldMem, class TimeStepperBuilder, LogicalToPhysicalMapping&gt;<br>_Owning manager for_ [_**ElementwisePhysicalAdvPseudoPhysicalFootFinder**_](classpolar__foot__finder__details_1_1ElementwisePhysicalAdvPseudoPhysicalFootFinder.md) _._ |



















































------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/advection/polar_foot_finders/physical_advection_pseudo_physical_foot_finder.hpp`

