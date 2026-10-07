

# File multipatch\_field.hpp



[**FileList**](files.md) **>** [**data\_types**](dir_2cbcac1ff0802c0a6551cceb4db325f2.md) **>** [**multipatch\_field.hpp**](multipatch__field_8hpp.md)

[Go to the source code of this file](multipatch__field_8hpp_source.md)



* `#include "multipatch_type.hpp"`













## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**ddcHelper**](namespaceddcHelper.md) <br> |




## Public Types

| Type | Name |
| ---: | :--- |
| typedef detail::MultipatchField&lt; ddc::detail::TypeSeq&lt; Patches... &gt;, ddc::detail::TypeSeq&lt; T&lt; Patches &gt;... &gt; &gt; | [**MultipatchField**](#typedef-multipatchfield)  <br>_A class to store field objects on patches._  |




## Public Attributes

| Type | Name |
| ---: | :--- |
|  constexpr bool | [**enable\_multipatch\_field**](#variable-enable_multipatch_field)   = `false`<br> |
|  constexpr bool | [**is\_multipatch\_field\_v**](#variable-is_multipatch_field_v)   = `enable\_multipatch\_field&lt;std::remove\_const\_t&lt;std::remove\_reference\_t&lt;T&gt;&gt;&gt;`<br> |












































## Public Types Documentation




### typedef MultipatchField 

_A class to store field objects on patches._ 
```C++
using MultipatchField =  detail:: MultipatchField<ddc::detail::TypeSeq<Patches...>, ddc::detail::TypeSeq<T<Patches>...> >;
```



See detail::MultipatchField for more details. This alias creates the type T&lt;Patch&gt; for each of the patches. Two type templates which return the same type for each patch therefore lead to the same MultipatchField.




**Template parameters:**


* `T` The type of the fields/derivative fields/vector fields that are stored on the given patches. 
* `Patches` The patches of the objects in the same order of the patches that the given objects are defined on. 




        

<hr>
## Public Attributes Documentation




### variable enable\_multipatch\_field 

```C++
constexpr bool enable_multipatch_field;
```




<hr>



### variable is\_multipatch\_field\_v 

```C++
constexpr bool is_multipatch_field_v;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/data_types/multipatch_field.hpp`

