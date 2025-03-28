

# File vector\_mapper.hpp



[**FileList**](files.md) **>** [**mapping**](dir_5300298560c4bf255ab9f36681603d89.md) **>** [**vector\_mapper.hpp**](vector__mapper_8hpp.md)

[Go to the source code of this file](vector__mapper_8hpp_source.md)



* `#include "vector_field.hpp"`
* `#include "vector_field_mem.hpp"`
* `#include "vector_index_tools.hpp"`
* `#include "view.hpp"`















## Classes

| Type | Name |
| ---: | :--- |
| class | [**VectorMapper&lt; VectorIndexSet&lt; XIn, YIn &gt;, VectorIndexSet&lt; XOut, YOut &gt;, Mapping, ExecSpace &gt;**](classVectorMapper_3_01VectorIndexSet_3_01XIn_00_01YIn_01_4_00_01VectorIndexSet_3_01XOut_00_01YOu77c12468788509067d2c0ef34f5e389c.md) &lt;class XIn, class YIn, class XOut, class YOut, class Mapping, class ExecSpace&gt;<br>_A class to map vector fields from one coordinate system to another._  |






















## Public Functions

| Type | Name |
| ---: | :--- |
|  auto | [**create\_geometry\_mirror\_view**](#function-create_geometry_mirror_view) (ExecSpace exec\_space, [**VectorField**](classVectorField.md)&lt; ElementType, IdxRangeType, VectorIndexSet&lt; [**X**](structX.md), [**Y**](structY.md) &gt;, typename ExecSpace::memory\_space, LayoutStridedPolicy &gt; vector\_field, Mapping mapping) <br>_A helper class to get a vector field on a pseudo-Cartesian geometry. If the pseudo-Cartesian geometry is the same as the Cartesian geometry then the same vector field is returned. If the pseudo-Cartesian geometry is different then the vectors in the vector field are mapped to the new geometry and a_ [_**VectorFieldMem**_](classVectorFieldMem.md) _is returned._ |




























## Public Functions Documentation




### function create\_geometry\_mirror\_view 

_A helper class to get a vector field on a pseudo-Cartesian geometry. If the pseudo-Cartesian geometry is the same as the Cartesian geometry then the same vector field is returned. If the pseudo-Cartesian geometry is different then the vectors in the vector field are mapped to the new geometry and a_ [_**VectorFieldMem**_](classVectorFieldMem.md) _is returned._
```C++
template<class ExecSpace, class Mapping, class ElementType, class IdxRangeType, class X, class Y, class LayoutStridedPolicy>
auto create_geometry_mirror_view (
    ExecSpace exec_space,
    VectorField < ElementType, IdxRangeType, VectorIndexSet< X , Y >, typename ExecSpace::memory_space, LayoutStridedPolicy > vector_field,
    Mapping mapping
) 
```





**Parameters:**


* `exec_space` The space on which the function is executed (CPU/GPU). 
* `vector_field` The vector field to be mapped to the pseudo-Cartesian geometry. 
* `mapping` A mapping describing the relation between the Cartesian and pseudo-Cartesian geometries.



**Returns:**

A [**VectorField**](classVectorField.md) or [**VectorFieldMem**](classVectorFieldMem.md) containing the vectors in the pseudo-Cartesian geometry. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/mapping/vector_mapper.hpp`

