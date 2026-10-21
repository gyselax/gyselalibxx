

# File metric\_tensor\_evaluator.hpp



[**FileList**](files.md) **>** [**coord\_transformations**](dir_67161c4ffadea73fddf46ea451c2f62c.md) **>** [**metric\_tensor\_evaluator.hpp**](metric__tensor__evaluator_8hpp.md)

[Go to the source code of this file](metric__tensor__evaluator_8hpp_source.md)



* `#include <concepts>`
* `#include <type_traits>`
* `#include "coord_transformation_tools.hpp"`
* `#include "ddc_aliases.hpp"`
* `#include "indexed_tensor.hpp"`
* `#include "inverse_jacobian_matrix.hpp"`
* `#include "tensor.hpp"`
* `#include "vector_index_tools.hpp"`
* `#include "view.hpp"`













## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**concepts**](namespaceconcepts.md) <br> |
| namespace | [**metric\_tensor\_detail**](namespacemetric__tensor__detail.md) <br> |


## Classes

| Type | Name |
| ---: | :--- |
| class | [**MetricTensorEvaluator**](classMetricTensorEvaluator.md) &lt;Mapping, class PositionCoordinate&gt;<br>_An operator for calculating the metric tensor._  |
| struct | [**IdxPlaceholder&lt; Coord&lt; Dims... &gt; &gt;**](structmetric__tensor__detail_1_1IdxPlaceholder_3_01Coord_3_01Dims_8_8_8_01_4_01_4.md) &lt;Dims&gt;<br>_A struct to obtain a Idx containing PlaceholderGrids associated with a coordinate._  |
| struct | [**PlaceholderGrid**](structmetric__tensor__detail_1_1PlaceholderGrid.md) &lt;class Dim&gt;<br>_A placeholder grid used to check that a metric tensor can be evaluated at an index. It is never instantiated in a discrete space._  |



















































------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/coord_transformations/metric_tensor_evaluator.hpp`

