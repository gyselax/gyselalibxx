

# File discrete\_mapping\_builder.hpp



[**FileList**](files.md) **>** [**mapping**](dir_5300298560c4bf255ab9f36681603d89.md) **>** [**discrete\_mapping\_builder.hpp**](discrete__mapping__builder_8hpp.md)

[Go to the source code of this file](discrete__mapping__builder_8hpp_source.md)



* `#include <ddc/ddc.hpp>`
* `#include <ddc/kernels/splines.hpp>`
* `#include "ddc_alias_inline_functions.hpp"`
* `#include "ddc_aliases.hpp"`
* `#include "discrete_to_cartesian.hpp"`















## Classes

| Type | Name |
| ---: | :--- |
| class | [**DiscreteToCartesianBuilder**](classDiscreteToCartesianBuilder.md) &lt;class [**X**](structX.md), class [**Y**](structY.md), class SplineBuilder, class SplineEvaluator&gt;<br>_A class to create a_ [_**DiscreteToCartesian**_](classDiscreteToCartesian.md) _instance from an analytical mapping. This class creates and stores splines memory spaces describing the analytical mapping. The discrete mapping is then created using the splines without copying data._ |
| class | [**RefinedDiscreteToCartesianBuilder**](classRefinedDiscreteToCartesianBuilder.md) &lt;class [**X**](structX.md), class [**Y**](structY.md), class SplineBuilder, class SplineEvaluator, ncells\_r, ncells\_theta&gt;<br>_A class to create a_ [_**DiscreteToCartesian**_](classDiscreteToCartesian.md) _instance from an analytical mapping. This class creates an instance which uses more refined splines than the provided builder and evaluator. This class creates and stores splines memory spaces describing the analytical mapping. The discrete mapping is then created using the splines without copying data._ |
| struct | [**BSplinesRRefined**](structRefinedDiscreteToCartesianBuilder_1_1BSplinesRRefined.md) <br>_The type of the radial B-splines on which the new mapping will be defined._  |
| struct | [**BSplinesThetaRefined**](structRefinedDiscreteToCartesianBuilder_1_1BSplinesThetaRefined.md) <br>_The type of the poloidal B-splines on which the new mapping will be defined._  |
| struct | [**GridRRefined**](structRefinedDiscreteToCartesianBuilder_1_1GridRRefined.md) <br>_The type of the grid of radial points on which the new mapping will be defined._  |
| struct | [**GridThetaRefined**](structRefinedDiscreteToCartesianBuilder_1_1GridThetaRefined.md) <br>_The type of the grid of poloidal points on which the new mapping will be defined._  |



















































------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/mapping/discrete_mapping_builder.hpp`

