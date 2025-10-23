

# File vector\_mapper.hpp



[**FileList**](files.md) **>** [**coord\_transformations**](dir_67161c4ffadea73fddf46ea451c2f62c.md) **>** [**vector\_mapper.hpp**](vector__mapper_8hpp.md)

[Go to the source code of this file](vector__mapper_8hpp_source.md)



* `#include "inverse_jacobian_matrix.hpp"`
* `#include "metric_tensor_evaluator.hpp"`
* `#include "vector_field.hpp"`
* `#include "vector_field_mem.hpp"`
* `#include "vector_index_tools.hpp"`
* `#include "view.hpp"`





































## Public Functions

| Type | Name |
| ---: | :--- |
|  void | [**copy\_to\_vector\_space**](#function-copy_to_vector_space) (ExecSpace exec\_space, [**VectorField**](classVectorField.md)&lt; ElementType, IdxRangeType, OutVectorSpace, typename ExecSpace::memory\_space, LayoutStridedPolicy &gt; vector\_field\_out, Mapping mapping, [**VectorConstField**](classVectorField.md)&lt; ElementType, IdxRangeType, InVectorSpace, typename ExecSpace::memory\_space, LayoutStridedPolicy &gt; vector\_field) <br>_A helper method to get a vector field on a different vector space._  |
|  auto | [**create\_mirror\_view\_and\_copy\_on\_vector\_space**](#function-create_mirror_view_and_copy_on_vector_space) (ExecSpace exec\_space, [**VectorField**](classVectorField.md)&lt; ElementType, IdxRangeType, InVectorSpace, typename ExecSpace::memory\_space, LayoutStridedPolicy &gt; vector\_field, Mapping mapping) <br>_A helper method to get a vector field on a different vector space. If the requested vector space is the same as the current vector space then the same vector field is returned. If the vector space is different then the vectors in the vector field are mapped to the new vector space and a_ [_**VectorFieldMem**_](classVectorFieldMem.md) _is returned._ |
|  KOKKOS\_INLINE\_FUNCTION [**Tensor**](classTensor.md)&lt; ElementType, OutVectorSpace &gt; | [**to\_vector\_space**](#function-to_vector_space) (Mapping const & mapping, CoordType const & coord, [**Tensor**](classTensor.md)&lt; ElementType, InVectorSpace &gt; const & in\_vector) <br>_A helper method to get a vector on a different vector space._  |




























## Public Functions Documentation




### function copy\_to\_vector\_space 

_A helper method to get a vector field on a different vector space._ 
```C++
template<class OutVectorSpace, class ExecSpace, class Mapping, class ElementType, class IdxRangeType, class InVectorSpace, class LayoutStridedPolicy>
void copy_to_vector_space (
    ExecSpace exec_space,
    VectorField < ElementType, IdxRangeType, OutVectorSpace, typename ExecSpace::memory_space, LayoutStridedPolicy > vector_field_out,
    Mapping mapping,
    VectorConstField < ElementType, IdxRangeType, InVectorSpace, typename ExecSpace::memory_space, LayoutStridedPolicy > vector_field
) 
```





**Parameters:**


* `exec_space` The space on which the function is executed (CPU/GPU). 
* `vector_field_out` The vector field where the vectors on the new vector space should be saved. 
* `mapping` A mapping describing the relation between the 2 vector spaces (or describing the relation between the vector space and the Cartesian space if a change of variance (covariant &lt;-&gt; contravariant) is required). 
* `vector_field` The vector field to be mapped to the new vector space.



**Returns:**

A [**VectorField**](classVectorField.md) or [**VectorFieldMem**](classVectorFieldMem.md) containing the vectors in the requested vector space. 





        

<hr>



### function create\_mirror\_view\_and\_copy\_on\_vector\_space 

_A helper method to get a vector field on a different vector space. If the requested vector space is the same as the current vector space then the same vector field is returned. If the vector space is different then the vectors in the vector field are mapped to the new vector space and a_ [_**VectorFieldMem**_](classVectorFieldMem.md) _is returned._
```C++
template<class OutVectorSpace, class ExecSpace, class Mapping, class ElementType, class IdxRangeType, class InVectorSpace, class LayoutStridedPolicy>
auto create_mirror_view_and_copy_on_vector_space (
    ExecSpace exec_space,
    VectorField < ElementType, IdxRangeType, InVectorSpace, typename ExecSpace::memory_space, LayoutStridedPolicy > vector_field,
    Mapping mapping
) 
```





**Parameters:**


* `exec_space` The space on which the function is executed (CPU/GPU). 
* `vector_field` The vector field to be mapped to the new vector space. 
* `mapping` A mapping describing the relation between the 2 vector spaces (or describing the relation between the vector space and the Cartesian space if a change of variance (covariant &lt;-&gt; contravariant) is required).



**Template parameters:**


* `OutVectorSpace` The vector space of the final vector.



**Returns:**

A [**VectorField**](classVectorField.md) or [**VectorFieldMem**](classVectorFieldMem.md) containing the vectors in the requested vector space. 





        

<hr>



### function to\_vector\_space 

_A helper method to get a vector on a different vector space._ 
```C++
template<class OutVectorSpace, class Mapping, class CoordType, class ElementType, class InVectorSpace>
KOKKOS_INLINE_FUNCTION Tensor < ElementType, OutVectorSpace > to_vector_space (
    Mapping const & mapping,
    CoordType const & coord,
    Tensor < ElementType, InVectorSpace > const & in_vector
) 
```





**Parameters:**


* `mapping` A mapping describing the relation between the 2 vector spaces (or describing the relation between the vector space and the Cartesian space if a change of variance (covariant &lt;-&gt; contravariant) is required). 
* `coord` The coordinate where the vector is defined. 
* `in_vector` The vector on the original vector space.



**Template parameters:**


* `OutVectorSpace` The vector space of the final vector.



**Returns:**

A [**VectorField**](classVectorField.md) or [**VectorFieldMem**](classVectorFieldMem.md) containing the vectors in the requested vector space. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/coord_transformations/vector_mapper.hpp`

