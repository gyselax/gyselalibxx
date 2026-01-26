

# File coord\_transformation\_tools.hpp



[**FileList**](files.md) **>** [**coord\_transformations**](dir_67161c4ffadea73fddf46ea451c2f62c.md) **>** [**coord\_transformation\_tools.hpp**](coord__transformation__tools_8hpp.md)

[Go to the source code of this file](coord__transformation__tools_8hpp_source.md)



* `#include <array>`
* `#include <concepts>`
* `#include <type_traits>`
* `#include <utility>`
* `#include "tensor.hpp"`
* `#include "vector_index_tools.hpp"`
* `#include "view.hpp"`













## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**concepts**](namespaceconcepts.md) <br> |




## Public Types

| Type | Name |
| ---: | :--- |
| typedef decltype(std::declval&lt; Mapping &gt;().get\_inverse\_mapping()) | [**inverse\_mapping\_t**](#typedef-inverse_mapping_t)  <br> |






## Public Static Attributes

| Type | Name |
| ---: | :--- |
|  constexpr bool | [**has\_singular\_o\_point\_inv\_jacobian\_v**](#variable-has_singular_o_point_inv_jacobian_v)   = `/* multi line expression */`<br> |
|  constexpr bool | [**is\_accessible\_v**](#variable-is_accessible_v)   = `/* multi line expression */`<br> |
|  constexpr bool | [**is\_coord\_transform\_with\_o\_point\_v**](#variable-is_coord_transform_with_o_point_v)   = `mapping\_detail::HasOPoint&lt;std::remove\_const\_t&lt;std::remove\_reference\_t&lt;Mapping&gt;&gt;&gt;::value`<br>_Indicates that a coordinate change operator is 2D with a curvilinear mapping showing an O-point._  |










































## Public Types Documentation




### typedef inverse\_mapping\_t 

```C++
using inverse_mapping_t =  decltype(std::declval<Mapping>().get_inverse_mapping());
```




<hr>
## Public Static Attributes Documentation




### variable has\_singular\_o\_point\_inv\_jacobian\_v 

```C++
constexpr bool has_singular_o_point_inv_jacobian_v;
```




<hr>



### variable is\_accessible\_v 

```C++
constexpr bool is_accessible_v;
```




<hr>



### variable is\_coord\_transform\_with\_o\_point\_v 

_Indicates that a coordinate change operator is 2D with a curvilinear mapping showing an O-point._ 
```C++
constexpr bool is_coord_transform_with_o_point_v;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/coord_transformations/coord_transformation_tools.hpp`

