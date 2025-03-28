

# File ddc\_alias\_inline\_functions.hpp



[**FileList**](files.md) **>** [**src**](dir_68267d1309a1af8e8297ef4c3efbcdba.md) **>** [**utils**](dir_313caf1132e152dd9b58bea13a4052ca.md) **>** [**ddc\_alias\_inline\_functions.hpp**](ddc__alias__inline__functions_8hpp.md)

[Go to the source code of this file](ddc__alias__inline__functions_8hpp_source.md)



* `#include <type_traits>`
* `#include <utility>`
* `#include <ddc/ddc.hpp>`
* `#include "ddc_aliases.hpp"`





















## Public Attributes

| Type | Name |
| ---: | :--- |
|  constexpr bool | [**enable\_data\_access\_methods**](#variable-enable_data_access_methods)   = `false`<br> |
|  constexpr bool | [**enable\_data\_access\_methods&lt; Field&lt; ElementType, IdxRangeType, LayoutType, MemoryType &gt; &gt;**](#variable-enable_data_access_methods-field-elementtype-idxrangetype-layouttype-memorytype)   = `true`<br> |
|  constexpr bool | [**enable\_data\_access\_methods&lt; FieldMem&lt; ElementType, IdxRangeType, MemoryType &gt; &gt;**](#variable-enable_data_access_methods-fieldmem-elementtype-idxrangetype-memorytype)   = `true`<br> |
|  constexpr bool | [**enable\_mem\_type**](#variable-enable_mem_type)   = `false`<br> |
|  constexpr bool | [**enable\_mem\_type&lt; FieldMem&lt; ElementType, IdxRangeType, MemoryType &gt; &gt;**](#variable-enable_mem_type-fieldmem-elementtype-idxrangetype-memorytype)   = `true`<br> |
|  constexpr bool | [**has\_data\_access\_methods\_v**](#variable-has_data_access_methods_v)   = `enable\_data\_access\_methods&lt;std::remove\_const\_t&lt;std::remove\_reference\_t&lt;Type&gt;&gt;&gt;`<br> |
|  constexpr bool | [**is\_mem\_type\_v**](#variable-is_mem_type_v)   = `enable\_mem\_type&lt;std::remove\_const\_t&lt;std::remove\_reference\_t&lt;Type&gt;&gt;&gt;`<br> |


## Public Static Attributes

| Type | Name |
| ---: | :--- |
|  constexpr bool | [**has\_idx\_range\_v**](#variable-has_idx_range_v)   = `detail::HasIdxRange&lt;Type&gt;::value`<br> |
|  constexpr bool | [**is\_gslx\_field\_v**](#variable-is_gslx_field_v)   = `detail::IsGslxField&lt;Type&gt;::value`<br> |














## Public Functions

| Type | Name |
| ---: | :--- |
|  KOKKOS\_INLINE\_FUNCTION auto | [**get\_const\_field**](#function-get_const_field) (FieldType && field) <br> |
|  auto | [**get\_const\_field**](#function-get_const_field) (FieldType && field) <br> |
|  KOKKOS\_INLINE\_FUNCTION auto | [**get\_field**](#function-get_field) (FieldType && field) <br> |
|  auto | [**get\_field**](#function-get_field) (FieldType && field) <br> |
|  auto | [**get\_idx\_range**](#function-get_idx_range) (FieldType const & field) noexcept<br> |
|  auto | [**get\_spline\_idx\_range**](#function-get_spline_idx_range) (SplineBuilder const & builder) noexcept<br> |




























## Public Attributes Documentation




### variable enable\_data\_access\_methods 

```C++
constexpr bool enable_data_access_methods;
```




<hr>



### variable enable\_data\_access\_methods&lt; Field&lt; ElementType, IdxRangeType, LayoutType, MemoryType &gt; &gt; 

```C++
constexpr bool enable_data_access_methods< Field< ElementType, IdxRangeType, LayoutType, MemoryType > >;
```




<hr>



### variable enable\_data\_access\_methods&lt; FieldMem&lt; ElementType, IdxRangeType, MemoryType &gt; &gt; 

```C++
constexpr bool enable_data_access_methods< FieldMem< ElementType, IdxRangeType, MemoryType > >;
```




<hr>



### variable enable\_mem\_type 

```C++
constexpr bool enable_mem_type;
```




<hr>



### variable enable\_mem\_type&lt; FieldMem&lt; ElementType, IdxRangeType, MemoryType &gt; &gt; 

```C++
constexpr bool enable_mem_type< FieldMem< ElementType, IdxRangeType, MemoryType > >;
```




<hr>



### variable has\_data\_access\_methods\_v 

```C++
constexpr bool has_data_access_methods_v;
```




<hr>



### variable is\_mem\_type\_v 

```C++
constexpr bool is_mem_type_v;
```




<hr>
## Public Static Attributes Documentation




### variable has\_idx\_range\_v 

```C++
constexpr bool has_idx_range_v;
```




<hr>



### variable is\_gslx\_field\_v 

```C++
constexpr bool is_gslx_field_v;
```




<hr>
## Public Functions Documentation




### function get\_const\_field 

```C++
template<class FieldType, std::enable_if_t<!is_mem_type_v< FieldType >, bool >>
KOKKOS_INLINE_FUNCTION auto get_const_field (
    FieldType && field
) 
```



A helper function to get a constant field from a FieldMem without allocating additional memory.




**Parameters:**


* `field` The field memory object. 



**Returns:**

The constant field. 





        

<hr>



### function get\_const\_field 

```C++
template<class FieldType, std::enable_if_t< is_mem_type_v< FieldType >, bool >>
inline auto get_const_field (
    FieldType && field
) 
```




<hr>



### function get\_field 

```C++
template<class FieldType, std::enable_if_t<!is_mem_type_v< FieldType >, bool >>
KOKKOS_INLINE_FUNCTION auto get_field (
    FieldType && field
) 
```



A helper function to get a modifiable field from a FieldMem without allocating additional memory.




**Parameters:**


* `field` The field memory object. 



**Returns:**

The modifiable field. 





        

<hr>



### function get\_field 

```C++
template<class FieldType, std::enable_if_t< is_mem_type_v< FieldType >, bool >>
inline auto get_field (
    FieldType && field
) 
```




<hr>



### function get\_idx\_range 

```C++
template<class... QueryGrids, class FieldType>
auto get_idx_range (
    FieldType const & field
) noexcept
```



A function to get the range of valid indices that can be used to index this field.




**Parameters:**


* `field` The field whose indices are of interest.



**Returns:**

The index range. 





        

<hr>



### function get\_spline\_idx\_range 

```C++
template<class SplineBuilder>
inline auto get_spline_idx_range (
    SplineBuilder const & builder
) noexcept
```



A function to get the range of valid indices that can be used to index a set of b-splines that are compatible with this spline builder.




**Parameters:**


* `builder` The spline builder.



**Returns:**

The index range. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/utils/ddc_alias_inline_functions.hpp`

