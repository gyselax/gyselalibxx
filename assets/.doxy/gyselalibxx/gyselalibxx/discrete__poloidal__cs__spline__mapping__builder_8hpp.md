

# File discrete\_poloidal\_cs\_spline\_mapping\_builder.hpp



[**FileList**](files.md) **>** [**coord\_transformations**](dir_67161c4ffadea73fddf46ea451c2f62c.md) **>** [**discrete\_poloidal\_cs\_spline\_mapping\_builder.hpp**](discrete__poloidal__cs__spline__mapping__builder_8hpp.md)

[Go to the source code of this file](discrete__poloidal__cs__spline__mapping__builder_8hpp_source.md)



* `#include <ddc/ddc.hpp>`
* `#include <ddc/kernels/splines.hpp>`
* `#include "ddc_alias_inline_functions.hpp"`
* `#include "ddc_aliases.hpp"`
* `#include "discrete_poloidal_cs_spline_mapping.hpp"`















## Classes

| Type | Name |
| ---: | :--- |
| class | [**DiscretePoloidalCSSplineMappingBuilder**](classDiscretePoloidalCSSplineMappingBuilder.md) &lt;class [**X**](structX.md), class [**Y**](structY.md), class SplineBuilder, class SplineEvaluator&gt;<br>_A class to create a_ [_**DiscretePoloidalCSSplineMapping**_](classDiscretePoloidalCSSplineMapping.md) _instance from an analytical mapping. This class creates and stores splines memory spaces describing the analytical mapping. The discrete mapping is then created using the splines without copying data._ |
| class | [**RefinedDiscretePoloidalCSSplineMappingBuilder**](classRefinedDiscretePoloidalCSSplineMappingBuilder.md) &lt;class [**X**](structX.md), class [**Y**](structY.md), class SplineBuilder, class SplineEvaluator, ncells\_r, ncells\_theta&gt;<br>_A class to create a_ [_**DiscretePoloidalCSSplineMapping**_](classDiscretePoloidalCSSplineMapping.md) _instance from an analytical mapping. This class creates an instance which uses more refined splines than the provided builder and evaluator. This class creates and stores splines memory spaces describing the analytical mapping. The discrete mapping is then created using the splines without copying data._ |
| struct | [**BSplinesRRefined**](structRefinedDiscretePoloidalCSSplineMappingBuilder_1_1BSplinesRRefined.md) <br>_The type of the radial B-splines on which the new mapping will be defined._  |
| struct | [**BSplinesThetaRefined**](structRefinedDiscretePoloidalCSSplineMappingBuilder_1_1BSplinesThetaRefined.md) <br>_The type of the poloidal B-splines on which the new mapping will be defined._  |
| struct | [**GridRRefined**](structRefinedDiscretePoloidalCSSplineMappingBuilder_1_1GridRRefined.md) <br>_The type of the grid of radial points on which the new mapping will be defined._  |
| struct | [**GridThetaRefined**](structRefinedDiscretePoloidalCSSplineMappingBuilder_1_1GridThetaRefined.md) <br>_The type of the grid of poloidal points on which the new mapping will be defined._  |



















































------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/coord_transformations/discrete_poloidal_cs_spline_mapping_builder.hpp`

