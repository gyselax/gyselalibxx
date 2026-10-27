

# Namespace vector\_mapper\_detail



[**Namespace List**](namespaces.md) **>** [**vector\_mapper\_detail**](namespacevector__mapper__detail.md)




















## Classes

| Type | Name |
| ---: | :--- |
| struct | [**EvaluationCoord**](structvector__mapper__detail_1_1EvaluationCoord.md) &lt;class MappingOrMetric&gt;<br>_The type of the coordinate at which a mapping or a metric tensor is evaluated._  |
| struct | [**EvaluationCoord&lt; Mapping &gt;**](structvector__mapper__detail_1_1EvaluationCoord_3_01Mapping_01_4.md) &lt;Mapping&gt;<br>_The type of the coordinate at which the Jacobian of a mapping is evaluated._  |
| struct | [**EvaluationCoord&lt; MetricTensorType &gt;**](structvector__mapper__detail_1_1EvaluationCoord_3_01MetricTensorType_01_4.md) &lt;MetricTensorType&gt;<br>_The type of the coordinate at which a metric tensor is evaluated._  |


## Public Types

| Type | Name |
| ---: | :--- |
| typedef typename [**EvaluationCoord**](structvector__mapper__detail_1_1EvaluationCoord.md)&lt; MappingOrMetric &gt;::type | [**evaluation\_coord\_t**](#typedef-evaluation_coord_t)  <br>_The type of the coordinate at which a mapping or a metric tensor is evaluated._  |






## Public Static Attributes

| Type | Name |
| ---: | :--- |
|  bool | [**are\_contravariant\_v**](#variable-are_contravariant_v)   = `/* multi line expression */`<br>_True if both vector spaces are contravariant._  |
|  bool | [**is\_arg\_to\_result\_v**](#variable-is_arg_to_result_v)   = `/* multi line expression */`<br>_True if the vector is mapped from the argument basis of the mapping to its result basis._  |
|  bool | [**is\_coordinate\_change\_only\_v**](#variable-is_coordinate_change_only_v)   = `/* multi line expression */`<br>_True if the vector spaces are associated with different coordinate systems and have the same variance._  |
|  bool | [**is\_result\_to\_arg\_v**](#variable-is_result_to_arg_v)   = `/* multi line expression */`<br>_True if the vector is mapped from the result basis of the mapping to its argument basis._  |
|  bool | [**is\_variance\_change\_only\_v**](#variable-is_variance_change_only_v)   = `/* multi line expression */`<br>_True if the vector spaces are the co/contra variant bases associated with the same coordinate system._  |










































## Public Types Documentation




### typedef evaluation\_coord\_t 

_The type of the coordinate at which a mapping or a metric tensor is evaluated._ 
```C++
using vector_mapper_detail::evaluation_coord_t =  typename EvaluationCoord<MappingOrMetric>::type;
```




<hr>
## Public Static Attributes Documentation




### variable are\_contravariant\_v 

_True if both vector spaces are contravariant._ 
```C++
bool vector_mapper_detail::are_contravariant_v;
```




<hr>



### variable is\_arg\_to\_result\_v 

_True if the vector is mapped from the argument basis of the mapping to its result basis._ 
```C++
bool vector_mapper_detail::is_arg_to_result_v;
```




<hr>



### variable is\_coordinate\_change\_only\_v 

_True if the vector spaces are associated with different coordinate systems and have the same variance._ 
```C++
bool vector_mapper_detail::is_coordinate_change_only_v;
```




<hr>



### variable is\_result\_to\_arg\_v 

_True if the vector is mapped from the result basis of the mapping to its argument basis._ 
```C++
bool vector_mapper_detail::is_result_to_arg_v;
```




<hr>



### variable is\_variance\_change\_only\_v 

_True if the vector spaces are the co/contra variant bases associated with the same coordinate system._ 
```C++
bool vector_mapper_detail::is_variance_change_only_v;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/coord_transformations/vector_mapper.hpp`

