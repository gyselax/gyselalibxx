

# File derivative\_field.hpp



[**FileList**](files.md) **>** [**data\_types**](dir_eaa769653453aaefd8cc10e98e9bb3eb.md) **>** [**derivative\_field.hpp**](derivative__field_8hpp.md)

[Go to the source code of this file](derivative__field_8hpp_source.md)



* `#include <ddc/ddc.hpp>`
* `#include "ddc_aliases.hpp"`
* `#include "derivative_field_common.hpp"`













## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**ddcHelper**](namespaceddcHelper.md) <br> |


## Classes

| Type | Name |
| ---: | :--- |
| class | [**DerivField&lt; ElementType, IdxRange&lt; DDims... &gt;, MemorySpace, LayoutStridedPolicy &gt;**](classDerivField_3_01ElementType_00_01IdxRange_3_01DDims_8_8_8_01_4_00_01MemorySpace_00_01LayoutStridedPolicy_01_4.md) &lt;class ElementType, DDims, class MemorySpace, class LayoutStridedPolicy&gt;<br>_A class which holds references to chunks of memory describing a field and its derivatives._  |


## Public Types

| Type | Name |
| ---: | :--- |
| typedef [**DerivField**](classDerivField.md)&lt; ElementType const, SupportType, MemorySpace, LayoutStridedPolicy &gt; | [**DerivConstField**](#typedef-derivconstfield)  <br> |




## Public Attributes

| Type | Name |
| ---: | :--- |
|  constexpr bool | [**enable\_borrowed\_deriv\_field&lt; DerivField&lt; ElementType, SupportType, MemorySpace, LayoutStridedPolicy &gt; &gt;**](#variable-enable_borrowed_deriv_field-derivfield-elementtype-supporttype-memoryspace-layoutstridedpolicy)   = `true`<br> |
|  constexpr bool | [**enable\_data\_access\_methods&lt; DerivField&lt; ElementType, SupportType, MemorySpace, LayoutStridedPolicy &gt; &gt;**](#variable-enable_data_access_methods-derivfield-elementtype-supporttype-memoryspace-layoutstridedpolicy)   = `true`<br> |
|  constexpr bool | [**enable\_deriv\_field&lt; DerivField&lt; ElementType, SupportType, MemorySpace, LayoutStridedPolicy &gt; &gt;**](#variable-enable_deriv_field-derivfield-elementtype-supporttype-memoryspace-layoutstridedpolicy)   = `true`<br> |












































## Public Types Documentation




### typedef DerivConstField 

```C++
using DerivConstField =  DerivField<ElementType const, SupportType, MemorySpace, LayoutStridedPolicy>;
```




<hr>
## Public Attributes Documentation




### variable enable\_borrowed\_deriv\_field&lt; DerivField&lt; ElementType, SupportType, MemorySpace, LayoutStridedPolicy &gt; &gt; 

```C++
constexpr bool enable_borrowed_deriv_field< DerivField< ElementType, SupportType, MemorySpace, LayoutStridedPolicy > >;
```




<hr>



### variable enable\_data\_access\_methods&lt; DerivField&lt; ElementType, SupportType, MemorySpace, LayoutStridedPolicy &gt; &gt; 

```C++
constexpr bool enable_data_access_methods< DerivField< ElementType, SupportType, MemorySpace, LayoutStridedPolicy > >;
```




<hr>



### variable enable\_deriv\_field&lt; DerivField&lt; ElementType, SupportType, MemorySpace, LayoutStridedPolicy &gt; &gt; 

```C++
constexpr bool enable_deriv_field< DerivField< ElementType, SupportType, MemorySpace, LayoutStridedPolicy > >;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/data_types/derivative_field.hpp`

