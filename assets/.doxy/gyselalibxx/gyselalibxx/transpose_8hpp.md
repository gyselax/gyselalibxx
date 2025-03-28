

# File transpose.hpp



[**FileList**](files.md) **>** [**src**](dir_68267d1309a1af8e8297ef4c3efbcdba.md) **>** [**utils**](dir_313caf1132e152dd9b58bea13a4052ca.md) **>** [**transpose.hpp**](transpose_8hpp.md)

[Go to the source code of this file](transpose_8hpp_source.md)



* `#include "ddc_alias_inline_functions.hpp"`
* `#include "ddc_aliases.hpp"`
* `#include "type_seq_tools.hpp"`













## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**ddcHelper**](namespaceddcHelper.md) <br> |
























## Public Functions

| Type | Name |
| ---: | :--- |
|  Field&lt; ElementType, IdxRangeIn, MemorySpace, LayoutStridedPolicyOut &gt; | [**transpose\_layout**](#function-transpose_layout) (ExecSpace const & execution\_space, Field&lt; ElementType, IdxRangeIn, MemorySpace, LayoutStridedPolicyOut &gt; transposed\_field, ConstField&lt; ElementType, IdxRangeOut, MemorySpace, LayoutStridedPolicyIn &gt; field\_to\_transpose) <br>_Copy data from a view in one layout into a span in a transposed layout._  |




























## Public Functions Documentation




### function transpose\_layout 

_Copy data from a view in one layout into a span in a transposed layout._ 
```C++
template<class ExecSpace, class ElementType, class IdxRangeOut, class MemorySpace, class IdxRangeIn, class LayoutStridedPolicyIn, class LayoutStridedPolicyOut>
Field< ElementType, IdxRangeIn, MemorySpace, LayoutStridedPolicyOut > transpose_layout (
    ExecSpace const & execution_space,
    Field< ElementType, IdxRangeIn, MemorySpace, LayoutStridedPolicyOut > transposed_field,
    ConstField< ElementType, IdxRangeOut, MemorySpace, LayoutStridedPolicyIn > field_to_transpose
) 
```



Layouts are described by DDC's DiscreteDomains and two layouts are considered to be a transposition of one another if both domains describe data on the same physical dimensions.




**Parameters:**


* `execution_space` The execution space (Host/Device) where the code will run. 
* `transposed_field` The span describing the data object which the data will be copied into. 
* `field_to_transpose` The constant span describing the data object where the original data is found.



**Returns:**

The transposed\_field describing the data object which the data was copied into. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/utils/transpose.hpp`

