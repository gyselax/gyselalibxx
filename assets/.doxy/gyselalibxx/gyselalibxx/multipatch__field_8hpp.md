

# File multipatch\_field.hpp



[**FileList**](files.md) **>** [**data\_types**](dir_2cbcac1ff0802c0a6551cceb4db325f2.md) **>** [**multipatch\_field.hpp**](multipatch__field_8hpp.md)

[Go to the source code of this file](multipatch__field_8hpp_source.md)



* `#include "multipatch_type.hpp"`













## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**ddcHelper**](namespaceddcHelper.md) <br> |


## Classes

| Type | Name |
| ---: | :--- |
| class | [**MultipatchField**](classMultipatchField.md) &lt;T, Patches&gt;<br>_A class to store field objects on patches._  |






## Public Attributes

| Type | Name |
| ---: | :--- |
|  constexpr bool | [**enable\_data\_access\_methods&lt; MultipatchField&lt; T, Patches... &gt; &gt;**](#variable-enable_data_access_methods-multipatchfield-t-patches)   = `true`<br> |
|  constexpr bool | [**enable\_multipatch\_field**](#variable-enable_multipatch_field)   = `false`<br> |
|  constexpr bool | [**enable\_multipatch\_field&lt; MultipatchField&lt; T, Patches... &gt; &gt;**](#variable-enable_multipatch_field-multipatchfield-t-patches)   = `true`<br> |
|  constexpr bool | [**enable\_multipatch\_type&lt; MultipatchField&lt; T, Patches... &gt; &gt;**](#variable-enable_multipatch_type-multipatchfield-t-patches)   = `true`<br> |
|  constexpr bool | [**is\_multipatch\_field\_v**](#variable-is_multipatch_field_v)   = `enable\_multipatch\_field&lt;std::remove\_const\_t&lt;std::remove\_reference\_t&lt;[**T**](structT.md)&gt;&gt;&gt;`<br> |












































## Public Attributes Documentation




### variable enable\_data\_access\_methods&lt; MultipatchField&lt; T, Patches... &gt; &gt; 

```C++
constexpr bool enable_data_access_methods< MultipatchField< T, Patches... > >;
```




<hr>



### variable enable\_multipatch\_field 

```C++
constexpr bool enable_multipatch_field;
```




<hr>



### variable enable\_multipatch\_field&lt; MultipatchField&lt; T, Patches... &gt; &gt; 

```C++
constexpr bool enable_multipatch_field< MultipatchField< T, Patches... > >;
```




<hr>



### variable enable\_multipatch\_type&lt; MultipatchField&lt; T, Patches... &gt; &gt; 

```C++
constexpr bool enable_multipatch_type< MultipatchField< T, Patches... > >;
```




<hr>



### variable is\_multipatch\_field\_v 

```C++
constexpr bool is_multipatch_field_v;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/data_types/multipatch_field.hpp`

