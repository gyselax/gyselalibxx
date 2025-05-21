

# File coord\_transformation\_tools.hpp



[**FileList**](files.md) **>** [**coord\_transformations**](dir_67161c4ffadea73fddf46ea451c2f62c.md) **>** [**coord\_transformation\_tools.hpp**](coord__transformation__tools_8hpp.md)

[Go to the source code of this file](coord__transformation__tools_8hpp_source.md)



* `#include <array>`
* `#include <type_traits>`
* `#include <utility>`
* `#include "tensor.hpp"`
* `#include "vector_index_tools.hpp"`
* `#include "view.hpp"`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef decltype(std::declval&lt; Mapping &gt;().get\_inverse\_mapping()) | [**inverse\_mapping\_t**](#typedef-inverse_mapping_t)  <br> |






## Public Static Attributes

| Type | Name |
| ---: | :--- |
|  constexpr bool | [**has\_inv\_jacobian\_v**](#variable-has_inv_jacobian_v)   = `mapping\_detail::DefinesInvJacobian&lt;Mapping, !RaiseError&gt;::value`<br> |
|  constexpr bool | [**has\_jacobian\_v**](#variable-has_jacobian_v)   = `mapping\_detail::DefinesJacobian&lt;Mapping, !RaiseError&gt;::value`<br> |
|  constexpr bool | [**has\_singular\_o\_point\_inv\_jacobian\_v**](#variable-has_singular_o_point_inv_jacobian_v)   = `/* multi line expression */`<br> |
|  constexpr bool | [**is\_accessible\_v**](#variable-is_accessible_v)   = `/* multi line expression */`<br> |
|  constexpr bool | [**is\_analytical\_mapping\_v**](#variable-is_analytical_mapping_v)   = `mapping\_detail::IsAnalyticalMapping&lt;Mapping&gt;::value`<br> |
|  constexpr bool | [**is\_curvilinear\_2d\_mapping\_v**](#variable-is_curvilinear_2d_mapping_v)   = `/* multi line expression */`<br> |
|  constexpr bool | [**is\_mapping\_v**](#variable-is_mapping_v)   = `mapping\_detail::IsMapping&lt;Mapping&gt;::value`<br> |










































## Public Types Documentation




### typedef inverse\_mapping\_t 

```C++
using inverse_mapping_t =  decltype(std::declval<Mapping>().get_inverse_mapping());
```




<hr>
## Public Static Attributes Documentation




### variable has\_inv\_jacobian\_v 

```C++
constexpr bool has_inv_jacobian_v;
```




<hr>



### variable has\_jacobian\_v 

```C++
constexpr bool has_jacobian_v;
```




<hr>



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



### variable is\_analytical\_mapping\_v 

```C++
constexpr bool is_analytical_mapping_v;
```




<hr>



### variable is\_curvilinear\_2d\_mapping\_v 

```C++
constexpr bool is_curvilinear_2d_mapping_v;
```




<hr>



### variable is\_mapping\_v 

```C++
constexpr bool is_mapping_v;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/coord_transformations/coord_transformation_tools.hpp`

