

# Namespace metric\_tensor\_detail



[**Namespace List**](namespaces.md) **>** [**metric\_tensor\_detail**](namespacemetric__tensor__detail.md)




















## Classes

| Type | Name |
| ---: | :--- |
| struct | [**IdxPlaceholder**](structmetric__tensor__detail_1_1IdxPlaceholder.md) &lt;class CoordType&gt;<br>_A struct to obtain a Idx containing PlaceholderGrids associated with a coordinate._  |
| struct | [**IdxPlaceholder&lt; Coord&lt; Dims... &gt; &gt;**](structmetric__tensor__detail_1_1IdxPlaceholder_3_01Coord_3_01Dims_8_8_8_01_4_01_4.md) &lt;Dims&gt;<br>_A struct to obtain a Idx containing PlaceholderGrids associated with a coordinate._  |
| struct | [**PlaceholderGrid**](structmetric__tensor__detail_1_1PlaceholderGrid.md) &lt;class Dim&gt;<br>_A placeholder grid used to check that a metric tensor can be evaluated at an index. It is never instantiated in a discrete space._  |


## Public Types

| Type | Name |
| ---: | :--- |
| typedef typename [**IdxPlaceholder**](structmetric__tensor__detail_1_1IdxPlaceholder.md)&lt; CoordType &gt;::type | [**placeholder\_idx\_t**](#typedef-placeholder_idx_t)  <br> |
















































## Public Types Documentation




### typedef placeholder\_idx\_t 

```C++
using metric_tensor_detail::placeholder_idx_t =  typename IdxPlaceholder<CoordType>::type;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/coord_transformations/metric_tensor_evaluator.hpp`

