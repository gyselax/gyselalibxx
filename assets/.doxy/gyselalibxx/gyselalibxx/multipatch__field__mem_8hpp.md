

# File multipatch\_field\_mem.hpp



[**FileList**](files.md) **>** [**data\_types**](dir_2cbcac1ff0802c0a6551cceb4db325f2.md) **>** [**multipatch\_field\_mem.hpp**](multipatch__field__mem_8hpp.md)

[Go to the source code of this file](multipatch__field__mem_8hpp_source.md)



* `#include "multipatch_field.hpp"`
* `#include "multipatch_type.hpp"`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef detail::MultipatchFieldMem&lt; ddc::detail::TypeSeq&lt; Patches... &gt;, ddc::detail::TypeSeq&lt; T&lt; Patches &gt;... &gt; &gt; | [**MultipatchFieldMem**](#typedef-multipatchfieldmem)  <br>_A class to store field memory block objects on patches._  |




## Public Attributes

| Type | Name |
| ---: | :--- |
|  constexpr bool | [**enable\_multipatch\_field\_mem**](#variable-enable_multipatch_field_mem)   = `false`<br> |
|  constexpr bool | [**is\_multipatch\_field\_mem\_v**](#variable-is_multipatch_field_mem_v)   = `enable\_multipatch\_field\_mem&lt;std::remove\_const\_t&lt;std::remove\_reference\_t&lt;T&gt;&gt;&gt;`<br> |












































## Public Types Documentation




### typedef MultipatchFieldMem 

_A class to store field memory block objects on patches._ 
```C++
using MultipatchFieldMem =  detail:: MultipatchFieldMem<ddc::detail::TypeSeq<Patches...>, ddc::detail::TypeSeq<T<Patches>...> >;
```



See detail::MultipatchFieldMem for more details. This alias creates the type T&lt;Patch&gt; for each of the patches. Two type templates which return the same type for each patch therefore lead to the same MultipatchFieldMem.




**Template parameters:**


* `T` The type of the FieldMem/DerivMem/VectorFieldMem that are stored on the given patches. 
* `Patches` The patches of the objects in the same order of the patches that the given objects are defined on. 




        

<hr>
## Public Attributes Documentation




### variable enable\_multipatch\_field\_mem 

```C++
constexpr bool enable_multipatch_field_mem;
```




<hr>



### variable is\_multipatch\_field\_mem\_v 

```C++
constexpr bool is_multipatch_field_mem_v;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/data_types/multipatch_field_mem.hpp`

