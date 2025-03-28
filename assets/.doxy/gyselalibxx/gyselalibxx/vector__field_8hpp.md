

# File vector\_field.hpp



[**FileList**](files.md) **>** [**data\_types**](dir_eaa769653453aaefd8cc10e98e9bb3eb.md) **>** [**vector\_field.hpp**](vector__field_8hpp.md)

[Go to the source code of this file](vector__field_8hpp_source.md)



* `#include "ddc_alias_inline_functions.hpp"`
* `#include "ddc_aliases.hpp"`
* `#include "ddc_helper.hpp"`
* `#include "vector_field_mem.hpp"`
* `#include "vector_index_tools.hpp"`













## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**ddcHelper**](namespaceddcHelper.md) <br> |


## Classes

| Type | Name |
| ---: | :--- |
| class | [**VectorField**](classVectorField.md) &lt;class ElementType, class IdxRangeType, class VectorIndexSetType, class MemorySpace, class LayoutStridedPolicy&gt;<br>_A class which holds multiple (scalar) fields in order to represent a vector field._  |


## Public Types

| Type | Name |
| ---: | :--- |
| typedef [**VectorConstField**](classVectorField.md)&lt; double, IdxRangeType, VectorIndexSetType, MemorySpace, LayoutStridedPolicy &gt; | [**DVectorConstField**](#typedef-dvectorconstfield)  <br> |
| typedef [**VectorField**](classVectorField.md)&lt; double, IdxRangeType, VectorIndexSetType, MemorySpace, LayoutStridedPolicy &gt; | [**DVectorField**](#typedef-dvectorfield)  <br> |
| typedef [**VectorField**](classVectorField.md)&lt; const ElementType, IdxRangeType, VectorIndexSetType, MemorySpace, LayoutStridedPolicy &gt; | [**VectorConstField**](#typedef-vectorconstfield)  <br> |




## Public Attributes

| Type | Name |
| ---: | :--- |
|  constexpr bool | [**enable\_borrowed\_vector\_field&lt; VectorField&lt; ElementType, IdxRangeType, VectorIndexSetType, MemorySpace, LayoutStridedPolicy &gt; &gt;**](#variable-enable_borrowed_vector_field-vectorfield-elementtype-idxrangetype-vectorindexsettype-memoryspace-layoutstridedpolicy)   = `true`<br> |
|  constexpr bool | [**enable\_data\_access\_methods&lt; VectorField&lt; ElementType, IdxRangeType, VectorIndexSetType, MemorySpace, LayoutStridedPolicy &gt; &gt;**](#variable-enable_data_access_methods-vectorfield-elementtype-idxrangetype-vectorindexsettype-memoryspace-layoutstridedpolicy)   = `true`<br> |
|  constexpr bool | [**enable\_vector\_field&lt; VectorField&lt; ElementType, IdxRangeType, VectorIndexSetType, MemorySpace, LayoutStridedPolicy &gt; &gt;**](#variable-enable_vector_field-vectorfield-elementtype-idxrangetype-vectorindexsettype-memoryspace-layoutstridedpolicy)   = `true`<br> |












































## Public Types Documentation




### typedef DVectorConstField 

```C++
using DVectorConstField =  VectorConstField< double, IdxRangeType, VectorIndexSetType, MemorySpace, LayoutStridedPolicy>;
```




<hr>



### typedef DVectorField 

```C++
using DVectorField =  VectorField<double, IdxRangeType, VectorIndexSetType, MemorySpace, LayoutStridedPolicy>;
```




<hr>



### typedef VectorConstField 

```C++
using VectorConstField =  VectorField< const ElementType, IdxRangeType, VectorIndexSetType, MemorySpace, LayoutStridedPolicy>;
```




<hr>
## Public Attributes Documentation




### variable enable\_borrowed\_vector\_field&lt; VectorField&lt; ElementType, IdxRangeType, VectorIndexSetType, MemorySpace, LayoutStridedPolicy &gt; &gt; 

```C++
constexpr bool enable_borrowed_vector_field< VectorField< ElementType, IdxRangeType, VectorIndexSetType, MemorySpace, LayoutStridedPolicy > >;
```




<hr>



### variable enable\_data\_access\_methods&lt; VectorField&lt; ElementType, IdxRangeType, VectorIndexSetType, MemorySpace, LayoutStridedPolicy &gt; &gt; 

```C++
constexpr bool enable_data_access_methods< VectorField< ElementType, IdxRangeType, VectorIndexSetType, MemorySpace, LayoutStridedPolicy > >;
```




<hr>



### variable enable\_vector\_field&lt; VectorField&lt; ElementType, IdxRangeType, VectorIndexSetType, MemorySpace, LayoutStridedPolicy &gt; &gt; 

```C++
constexpr bool enable_vector_field< VectorField< ElementType, IdxRangeType, VectorIndexSetType, MemorySpace, LayoutStridedPolicy > >;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/data_types/vector_field.hpp`

