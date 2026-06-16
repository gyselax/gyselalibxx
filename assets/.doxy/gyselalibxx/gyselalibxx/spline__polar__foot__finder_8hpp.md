

# File spline\_polar\_foot\_finder.hpp



[**FileList**](files.md) **>** [**advection**](dir_b90fde0f10c67a9aef841a6e6700f1f6.md) **>** [**spline\_polar\_foot\_finder.hpp**](spline__polar__foot__finder_8hpp.md)

[Go to the source code of this file](spline__polar__foot__finder_8hpp_source.md)



* `#include <functional>`
* `#include "circular_to_cartesian.hpp"`
* `#include "combined_mapping.hpp"`
* `#include "coord_transformation_tools.hpp"`
* `#include "ddc_alias_inline_functions.hpp"`
* `#include "ddc_aliases.hpp"`
* `#include "geometry_pseudo_cartesian.hpp"`
* `#include "itimestepper.hpp"`
* `#include "l_norm_tools.hpp"`
* `#include "vector_index_tools.hpp"`
* `#include "vector_mapper.hpp"`















## Classes

| Type | Name |
| ---: | :--- |
| class | [**ElementwiseSplinePolarFootFinder**](classElementwiseSplinePolarFootFinder.md) &lt;class [**GridR**](structGridR.md), class [**GridTheta**](structGridTheta.md), class IdxRangeOperator, class SplineRThetaEvaluatorAdvection, class PseudoPhysicalToAdvectionMapping, class PseudoPhysicalToLogicalMapping, class LogicalToPseudoPhysicalMapping, class AdvecCoefField, class TimeStepper&gt;<br>_A device-callable functor that finds the foot of the characteristic at a single grid point in polar coordinates using spline interpolation._  |
| class | [**ElementwiseSplinePolarFootFinderMem**](classElementwiseSplinePolarFootFinderMem.md) &lt;class [**GridR**](structGridR.md), class [**GridTheta**](structGridTheta.md), class IdxRangeOperator, class SplineRThetaEvaluatorAdvection, class PseudoPhysicalToAdvectionMapping, class PseudoPhysicalToLogicalMapping, class LogicalToPseudoPhysicalMapping, class AdvecCoefFieldMem, class TimeStepper&gt;<br>_The owning counterpart to_ [_**ElementwiseSplinePolarFootFinder**_](classElementwiseSplinePolarFootFinder.md) _._ |
| class | [**SplinePolarFootFinder**](classSplinePolarFootFinder.md) &lt;class IdxRangeBatched, class TimeStepperBuilder, LogicalToPhysicalMapping, LogicalToPseudoPhysicalMapping, class SplineRThetaBuilderAdvection, class SplineRThetaEvaluatorAdvection&gt;<br>_A class to find the foot of the characteristics on the_ \((r,\theta)\) _plane._ |



















































------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/advection/spline_polar_foot_finder.hpp`

