

# File vector\_mapper.hpp



[**FileList**](files.md) **>** [**coord\_transformations**](dir_67161c4ffadea73fddf46ea451c2f62c.md) **>** [**vector\_mapper.hpp**](vector__mapper_8hpp.md)

[Go to the source code of this file](vector__mapper_8hpp_source.md)



* `#include "inverse_jacobian_matrix.hpp"`
* `#include "metric_tensor_evaluator.hpp"`
* `#include "vector_field.hpp"`
* `#include "vector_field_mem.hpp"`
* `#include "vector_index_tools.hpp"`
* `#include "view.hpp"`













## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**vector\_mapper\_detail**](namespacevector__mapper__detail.md) <br> |


## Classes

| Type | Name |
| ---: | :--- |
| struct | [**EvaluationCoord&lt; Mapping &gt;**](structvector__mapper__detail_1_1EvaluationCoord_3_01Mapping_01_4.md) &lt;Mapping&gt;<br>_The type of the coordinate at which the Jacobian of a mapping is evaluated._  |
| struct | [**EvaluationCoord&lt; MetricTensorType &gt;**](structvector__mapper__detail_1_1EvaluationCoord_3_01MetricTensorType_01_4.md) &lt;MetricTensorType&gt;<br>_The type of the coordinate at which a metric tensor is evaluated._  |






















## Public Functions

| Type | Name |
| ---: | :--- |
|  void | [**copy\_to\_vector\_space**](#function-copy_to_vector_space) (ExecSpace exec\_space, [**VectorField**](classVectorField.md)&lt; ElementType, IdxRangeType, OutVectorSpace, typename ExecSpace::memory\_space, LayoutStridedPolicy &gt; vector\_field\_out, MappingOrMetric mapping\_or\_metric, [**VectorConstField**](classVectorField.md)&lt; ElementType, IdxRangeType, InVectorSpace, typename ExecSpace::memory\_space, LayoutStridedPolicy &gt; vector\_field) <br>_A helper method to get a vector field on a different vector space._  |
|  auto | [**create\_mirror\_view\_and\_copy\_on\_vector\_space**](#function-create_mirror_view_and_copy_on_vector_space) (std::string const & label, ExecSpace exec\_space, [**VectorField**](classVectorField.md)&lt; ElementType, IdxRangeType, InVectorSpace, typename ExecSpace::memory\_space, LayoutStridedPolicy &gt; vector\_field, MappingOrMetric mapping\_or\_metric) <br>_A helper method to get a vector field on a different vector space. If the requested vector space is the same as the current vector space then the same vector field is returned. If the vector space is different then the vectors in the vector field are mapped to the new vector space and a_ [_**VectorFieldMem**_](classVectorFieldMem.md) _is returned._ |
|  auto | [**create\_mirror\_view\_and\_copy\_on\_vector\_space**](#function-create_mirror_view_and_copy_on_vector_space) (ExecSpace exec\_space, [**VectorField**](classVectorField.md)&lt; ElementType, IdxRangeType, InVectorSpace, typename ExecSpace::memory\_space, LayoutStridedPolicy &gt; vector\_field, MappingOrMetric mapping\_or\_metric) <br>_A helper method to get a vector field on a different vector space. If the requested vector space is the same as the current vector space then the same vector field is returned. If the vector space is different then the vectors in the vector field are mapped to the new vector space and a_ [_**VectorFieldMem**_](classVectorFieldMem.md) _is returned._ |
|  KOKKOS\_INLINE\_FUNCTION [**Tensor**](classTensor.md)&lt; ElementType, OutVectorSpace &gt; | [**to\_vector\_space**](#function-to_vector_space) (MappingOrMetric const & mapping\_or\_metric, CoordType const & coord, [**Tensor**](classTensor.md)&lt; ElementType, InVectorSpace &gt; const & in\_vector) <br>_A helper method to get a vector on a different vector space. This overload handles the case where the vector is already defined on the requested vector space. It can be called with either a mapping or a metric tensor so that it is available to both the coordinate change and the variance change functions._  |
|  KOKKOS\_INLINE\_FUNCTION [**Tensor**](classTensor.md)&lt; ElementType, OutVectorSpace &gt; | [**to\_vector\_space**](#function-to_vector_space) (MetricTensorType const & metric\_tensor, CoordType const & coord, [**Tensor**](classTensor.md)&lt; ElementType, InVectorSpace &gt; const & in\_vector) <br>_A helper method to get a vector on a different vector space. This overload handles the case where a contravariant vector is mapped to the covariant basis associated with the same coordinate system._  |
|  KOKKOS\_INLINE\_FUNCTION [**Tensor**](classTensor.md)&lt; ElementType, OutVectorSpace &gt; | [**to\_vector\_space**](#function-to_vector_space) (MetricTensorType const & metric\_tensor, CoordType const & coord, [**Tensor**](classTensor.md)&lt; ElementType, InVectorSpace &gt; const & in\_vector) <br>_A helper method to get a vector on a different vector space. This overload handles the case where a covariant vector is mapped to the contravariant basis associated with the same coordinate system._  |
|  KOKKOS\_INLINE\_FUNCTION [**Tensor**](classTensor.md)&lt; ElementType, OutVectorSpace &gt; | [**to\_vector\_space**](#function-to_vector_space) (Mapping const & mapping, CoordType const & coord, [**Tensor**](classTensor.md)&lt; ElementType, InVectorSpace &gt; const & in\_vector) <br>_A helper method to get a vector on a different vector space. This overload handles the case where a contravariant vector is mapped from the argument basis of the mapping to its result basis._  |
|  KOKKOS\_INLINE\_FUNCTION [**Tensor**](classTensor.md)&lt; ElementType, OutVectorSpace &gt; | [**to\_vector\_space**](#function-to_vector_space) (Mapping const & mapping, CoordType const & coord, [**Tensor**](classTensor.md)&lt; ElementType, InVectorSpace &gt; const & in\_vector) <br>_A helper method to get a vector on a different vector space. This overload handles the case where a covariant vector is mapped from the argument basis of the mapping to its result basis._  |
|  KOKKOS\_INLINE\_FUNCTION [**Tensor**](classTensor.md)&lt; ElementType, OutVectorSpace &gt; | [**to\_vector\_space**](#function-to_vector_space) (Mapping const & mapping, CoordType const & coord, [**Tensor**](classTensor.md)&lt; ElementType, InVectorSpace &gt; const & in\_vector) <br>_A helper method to get a vector on a different vector space. This overload handles the case where a contravariant vector is mapped from the result basis of the mapping to its argument basis._  |
|  KOKKOS\_INLINE\_FUNCTION [**Tensor**](classTensor.md)&lt; ElementType, OutVectorSpace &gt; | [**to\_vector\_space**](#function-to_vector_space) (Mapping const & mapping, CoordType const & coord, [**Tensor**](classTensor.md)&lt; ElementType, InVectorSpace &gt; const & in\_vector) <br>_A helper method to get a vector on a different vector space. This overload handles the case where a covariant vector is mapped from the result basis of the mapping to its argument basis._  |




























## Public Functions Documentation




### function copy\_to\_vector\_space 

_A helper method to get a vector field on a different vector space._ 
```C++
template<class OutVectorSpace, class ExecSpace, class MappingOrMetric, class ElementType, class IdxRangeType, class InVectorSpace, class LayoutStridedPolicy>
void copy_to_vector_space (
    ExecSpace exec_space,
    VectorField < ElementType, IdxRangeType, OutVectorSpace, typename ExecSpace::memory_space, LayoutStridedPolicy > vector_field_out,
    MappingOrMetric mapping_or_metric,
    VectorConstField < ElementType, IdxRangeType, InVectorSpace, typename ExecSpace::memory_space, LayoutStridedPolicy > vector_field
) 
```





**Parameters:**


* `exec_space` The space on which the function is executed (CPU/GPU). 
* `vector_field_out` The vector field where the vectors on the new vector space should be saved. 
* `mapping_or_metric` A mapping describing the relation between the 2 vector spaces if the coordinate system changes, or a metric tensor if the variance (covariant &lt;-&gt; contravariant) changes. 
* `vector_field` The vector field to be mapped to the new vector space.



**Returns:**

A [**VectorField**](classVectorField.md) or [**VectorFieldMem**](classVectorFieldMem.md) containing the vectors in the requested vector space. 





        

<hr>



### function create\_mirror\_view\_and\_copy\_on\_vector\_space 

_A helper method to get a vector field on a different vector space. If the requested vector space is the same as the current vector space then the same vector field is returned. If the vector space is different then the vectors in the vector field are mapped to the new vector space and a_ [_**VectorFieldMem**_](classVectorFieldMem.md) _is returned._
```C++
template<class OutVectorSpace, class ExecSpace, class MappingOrMetric, class ElementType, class IdxRangeType, class InVectorSpace, class LayoutStridedPolicy>
auto create_mirror_view_and_copy_on_vector_space (
    std::string const & label,
    ExecSpace exec_space,
    VectorField < ElementType, IdxRangeType, InVectorSpace, typename ExecSpace::memory_space, LayoutStridedPolicy > vector_field,
    MappingOrMetric mapping_or_metric
) 
```





**Parameters:**


* `label` Label to use in case a new memory allocation. 
* `exec_space` The space on which the function is executed (CPU/GPU). 
* `vector_field` The vector field to be mapped to the new vector space. 
* `mapping_or_metric` A mapping describing the relation between the 2 vector spaces if the coordinate system changes, or a metric tensor if the variance (covariant &lt;-&gt; contravariant) changes.



**Template parameters:**


* `OutVectorSpace` The vector space of the final vector.



**Returns:**

A [**VectorField**](classVectorField.md) or [**VectorFieldMem**](classVectorFieldMem.md) containing the vectors in the requested vector space. 





        

<hr>



### function create\_mirror\_view\_and\_copy\_on\_vector\_space 

_A helper method to get a vector field on a different vector space. If the requested vector space is the same as the current vector space then the same vector field is returned. If the vector space is different then the vectors in the vector field are mapped to the new vector space and a_ [_**VectorFieldMem**_](classVectorFieldMem.md) _is returned._
```C++
template<class OutVectorSpace, class ExecSpace, class MappingOrMetric, class ElementType, class IdxRangeType, class InVectorSpace, class LayoutStridedPolicy>
auto create_mirror_view_and_copy_on_vector_space (
    ExecSpace exec_space,
    VectorField < ElementType, IdxRangeType, InVectorSpace, typename ExecSpace::memory_space, LayoutStridedPolicy > vector_field,
    MappingOrMetric mapping_or_metric
) 
```





**Parameters:**


* `exec_space` The space on which the function is executed (CPU/GPU). 
* `vector_field` The vector field to be mapped to the new vector space. 
* `mapping_or_metric` A mapping describing the relation between the 2 vector spaces if the coordinate system changes, or a metric tensor if the variance (covariant &lt;-&gt; contravariant) changes.



**Template parameters:**


* `OutVectorSpace` The vector space of the final vector.



**Returns:**

A [**VectorField**](classVectorField.md) or [**VectorFieldMem**](classVectorFieldMem.md) containing the vectors in the requested vector space. 





        

<hr>



### function to\_vector\_space 

_A helper method to get a vector on a different vector space. This overload handles the case where the vector is already defined on the requested vector space. It can be called with either a mapping or a metric tensor so that it is available to both the coordinate change and the variance change functions._ 
```C++
template<class OutVectorSpace, class MappingOrMetric, class CoordType, class ElementType, class InVectorSpace>
KOKKOS_INLINE_FUNCTION Tensor < ElementType, OutVectorSpace > to_vector_space (
    MappingOrMetric const & mapping_or_metric,
    CoordType const & coord,
    Tensor < ElementType, InVectorSpace > const & in_vector
) 
```





**Parameters:**


* `mapping_or_metric` A mapping or a metric tensor (unused). 
* `coord` The coordinate (or the index of the grid point) where the vector is defined. 
* `in_vector` The vector on the original vector space.



**Template parameters:**


* `OutVectorSpace` The vector space of the final vector.



**Returns:**

The vector on the requested vector space. 





        

<hr>



### function to\_vector\_space 

_A helper method to get a vector on a different vector space. This overload handles the case where a contravariant vector is mapped to the covariant basis associated with the same coordinate system._ 
```C++
template<class OutVectorSpace, concepts::MetricTensor MetricTensorType, class CoordType, class ElementType, class InVectorSpace>
KOKKOS_INLINE_FUNCTION Tensor < ElementType, OutVectorSpace > to_vector_space (
    MetricTensorType const & metric_tensor,
    CoordType const & coord,
    Tensor < ElementType, InVectorSpace > const & in_vector
) 
```





**Parameters:**


* `metric_tensor` An object which provides the metric tensor. 
* `coord` The coordinate (or the index of the grid point) where the vector is defined. 
* `in_vector` The vector on the original vector space.



**Template parameters:**


* `OutVectorSpace` The vector space of the final vector.



**Returns:**

The vector on the requested vector space. 





        

<hr>



### function to\_vector\_space 

_A helper method to get a vector on a different vector space. This overload handles the case where a covariant vector is mapped to the contravariant basis associated with the same coordinate system._ 
```C++
template<class OutVectorSpace, concepts::MetricTensor MetricTensorType, class CoordType, class ElementType, class InVectorSpace>
KOKKOS_INLINE_FUNCTION Tensor < ElementType, OutVectorSpace > to_vector_space (
    MetricTensorType const & metric_tensor,
    CoordType const & coord,
    Tensor < ElementType, InVectorSpace > const & in_vector
) 
```





**Parameters:**


* `metric_tensor` An object which provides the metric tensor. 
* `coord` The coordinate (or the index of the grid point) where the vector is defined. 
* `in_vector` The vector on the original vector space.



**Template parameters:**


* `OutVectorSpace` The vector space of the final vector.



**Returns:**

The vector on the requested vector space. 





        

<hr>



### function to\_vector\_space 

_A helper method to get a vector on a different vector space. This overload handles the case where a contravariant vector is mapped from the argument basis of the mapping to its result basis._ 
```C++
template<class OutVectorSpace, concepts::MappingWithJacobian Mapping, class CoordType, class ElementType, class InVectorSpace>
KOKKOS_INLINE_FUNCTION Tensor < ElementType, OutVectorSpace > to_vector_space (
    Mapping const & mapping,
    CoordType const & coord,
    Tensor < ElementType, InVectorSpace > const & in_vector
) 
```





**Parameters:**


* `mapping` A mapping describing the relation between the 2 vector spaces. 
* `coord` The coordinate where the vector is defined. 
* `in_vector` The vector on the original vector space.



**Template parameters:**


* `OutVectorSpace` The vector space of the final vector.



**Returns:**

The vector on the requested vector space. 





        

<hr>



### function to\_vector\_space 

_A helper method to get a vector on a different vector space. This overload handles the case where a covariant vector is mapped from the argument basis of the mapping to its result basis._ 
```C++
template<class OutVectorSpace, concepts::MappingWithJacobian Mapping, class CoordType, class ElementType, class InVectorSpace>
KOKKOS_INLINE_FUNCTION Tensor < ElementType, OutVectorSpace > to_vector_space (
    Mapping const & mapping,
    CoordType const & coord,
    Tensor < ElementType, InVectorSpace > const & in_vector
) 
```





**Parameters:**


* `mapping` A mapping describing the relation between the 2 vector spaces. 
* `coord` The coordinate where the vector is defined. 
* `in_vector` The vector on the original vector space.



**Template parameters:**


* `OutVectorSpace` The vector space of the final vector.



**Returns:**

The vector on the requested vector space. 





        

<hr>



### function to\_vector\_space 

_A helper method to get a vector on a different vector space. This overload handles the case where a contravariant vector is mapped from the result basis of the mapping to its argument basis._ 
```C++
template<class OutVectorSpace, concepts::MappingWithJacobian Mapping, class CoordType, class ElementType, class InVectorSpace>
KOKKOS_INLINE_FUNCTION Tensor < ElementType, OutVectorSpace > to_vector_space (
    Mapping const & mapping,
    CoordType const & coord,
    Tensor < ElementType, InVectorSpace > const & in_vector
) 
```





**Parameters:**


* `mapping` A mapping describing the relation between the 2 vector spaces. 
* `coord` The coordinate where the vector is defined. 
* `in_vector` The vector on the original vector space.



**Template parameters:**


* `OutVectorSpace` The vector space of the final vector.



**Returns:**

The vector on the requested vector space. 





        

<hr>



### function to\_vector\_space 

_A helper method to get a vector on a different vector space. This overload handles the case where a covariant vector is mapped from the result basis of the mapping to its argument basis._ 
```C++
template<class OutVectorSpace, concepts::MappingWithJacobian Mapping, class CoordType, class ElementType, class InVectorSpace>
KOKKOS_INLINE_FUNCTION Tensor < ElementType, OutVectorSpace > to_vector_space (
    Mapping const & mapping,
    CoordType const & coord,
    Tensor < ElementType, InVectorSpace > const & in_vector
) 
```





**Parameters:**


* `mapping` A mapping describing the relation between the 2 vector spaces. 
* `coord` The coordinate where the vector is defined. 
* `in_vector` The vector on the original vector space.



**Template parameters:**


* `OutVectorSpace` The vector space of the final vector.



**Returns:**

The vector on the requested vector space. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/coord_transformations/vector_mapper.hpp`

