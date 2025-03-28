

# File vector\_field\_mem.hpp



[**FileList**](files.md) **>** [**data\_types**](dir_eaa769653453aaefd8cc10e98e9bb3eb.md) **>** [**vector\_field\_mem.hpp**](vector__field__mem_8hpp.md)

[Go to the source code of this file](vector__field__mem_8hpp_source.md)



* `#include <ddc/ddc.hpp>`
* `#include "ddc_alias_inline_functions.hpp"`
* `#include "ddc_aliases.hpp"`
* `#include "ddc_helper.hpp"`
* `#include "vector_field_common.hpp"`
* `#include "vector_index_tools.hpp"`















## Classes

| Type | Name |
| ---: | :--- |
| class | [**VectorFieldMem**](classVectorFieldMem.md) &lt;class ElementType, class IdxRangeType, class VectorIndexSetType, class MemSpace&gt;<br>_Pre-declaration of_ [_**VectorFieldMem**_](classVectorFieldMem.md) _._ |


## Public Types

| Type | Name |
| ---: | :--- |
| typedef [**VectorFieldMem**](classVectorFieldMem.md)&lt; double, IdxRangeType, VectorIndexSetType, MemSpace &gt; | [**DVectorFieldMem**](#typedef-dvectorfieldmem)  <br> |




## Public Attributes

| Type | Name |
| ---: | :--- |
|  constexpr bool | [**enable\_data\_access\_methods&lt; VectorFieldMem&lt; ElementType, IdxRangeType, DimSeq, Allocator &gt; &gt;**](#variable-enable_data_access_methods-vectorfieldmem-elementtype-idxrangetype-dimseq-allocator)   = `true`<br> |
|  constexpr bool | [**enable\_mem\_type&lt; VectorFieldMem&lt; ElementType, IdxRangeType, DimSeq, Allocator &gt; &gt;**](#variable-enable_mem_type-vectorfieldmem-elementtype-idxrangetype-dimseq-allocator)   = `true`<br> |
|  constexpr bool | [**enable\_vector\_field&lt; VectorFieldMem&lt; ElementType, IdxRangeType, DimSeq, MemSpace &gt; &gt;**](#variable-enable_vector_field-vectorfieldmem-elementtype-idxrangetype-dimseq-memspace)   = `true`<br> |












































## Public Types Documentation




### typedef DVectorFieldMem 

```C++
using DVectorFieldMem =  VectorFieldMem<double, IdxRangeType, VectorIndexSetType, MemSpace>;
```




<hr>
## Public Attributes Documentation




### variable enable\_data\_access\_methods&lt; VectorFieldMem&lt; ElementType, IdxRangeType, DimSeq, Allocator &gt; &gt; 

```C++
constexpr bool enable_data_access_methods< VectorFieldMem< ElementType, IdxRangeType, DimSeq, Allocator > >;
```




<hr>



### variable enable\_mem\_type&lt; VectorFieldMem&lt; ElementType, IdxRangeType, DimSeq, Allocator &gt; &gt; 

```C++
constexpr bool enable_mem_type< VectorFieldMem< ElementType, IdxRangeType, DimSeq, Allocator > >;
```




<hr>



### variable enable\_vector\_field&lt; VectorFieldMem&lt; ElementType, IdxRangeType, DimSeq, MemSpace &gt; &gt; 

```C++
constexpr bool enable_vector_field< VectorFieldMem< ElementType, IdxRangeType, DimSeq, MemSpace > >;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/data_types/vector_field_mem.hpp`

