

# File mapping\_tools.hpp



[**FileList**](files.md) **>** [**mapping**](dir_5300298560c4bf255ab9f36681603d89.md) **>** [**mapping\_tools.hpp**](mapping__tools_8hpp.md)

[Go to the source code of this file](mapping__tools_8hpp_source.md)



* `#include <array>`
* `#include <type_traits>`
* `#include <utility>`
* `#include "view.hpp"`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef decltype(std::declval&lt; Mapping &gt;().get\_inverse\_mapping()) | [**inverse\_mapping\_t**](#typedef-inverse_mapping_t)  <br> |






## Public Static Attributes

| Type | Name |
| ---: | :--- |
|  constexpr bool | [**has\_2d\_inv\_jacobian\_v**](#variable-has_2d_inv_jacobian_v)   = `mapping\_detail::Defines2DInvJacobian&lt;Mapping, CoordinateType&gt;::value`<br> |
|  constexpr bool | [**has\_2d\_jacobian\_v**](#variable-has_2d_jacobian_v)   = `mapping\_detail::Defines2DJacobian&lt;Mapping, CoordinateType&gt;::value`<br> |
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




### variable has\_2d\_inv\_jacobian\_v 

```C++
constexpr bool has_2d_inv_jacobian_v;
```




<hr>



### variable has\_2d\_jacobian\_v 

```C++
constexpr bool has_2d_jacobian_v;
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
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/mapping/mapping_tools.hpp`

