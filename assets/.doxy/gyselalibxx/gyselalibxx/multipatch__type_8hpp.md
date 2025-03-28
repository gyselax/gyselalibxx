

# File multipatch\_type.hpp



[**FileList**](files.md) **>** [**data\_types**](dir_2cbcac1ff0802c0a6551cceb4db325f2.md) **>** [**multipatch\_type.hpp**](multipatch__type_8hpp.md)

[Go to the source code of this file](multipatch__type_8hpp_source.md)



* `#include "ddc_alias_inline_functions.hpp"`
* `#include "ddc_aliases.hpp"`
* `#include "types.hpp"`















## Classes

| Type | Name |
| ---: | :--- |
| class | [**MultipatchType**](classMultipatchType.md) &lt;T, Patches&gt;<br>_A class to store several objects that are of a type which is templated by the patch._  |






## Public Attributes

| Type | Name |
| ---: | :--- |
|  constexpr bool | [**enable\_multipatch\_type**](#variable-enable_multipatch_type)   = `false`<br> |
|  constexpr bool | [**enable\_multipatch\_type&lt; MultipatchType&lt; T, Patches... &gt; &gt;**](#variable-enable_multipatch_type-multipatchtype-t-patches)   = `true`<br> |
|  constexpr bool | [**is\_multipatch\_type\_v**](#variable-is_multipatch_type_v)   = `enable\_multipatch\_type&lt;std::remove\_const\_t&lt;std::remove\_reference\_t&lt;[**T**](structT.md)&gt;&gt;&gt;`<br> |












































## Public Attributes Documentation




### variable enable\_multipatch\_type 

```C++
constexpr bool enable_multipatch_type;
```




<hr>



### variable enable\_multipatch\_type&lt; MultipatchType&lt; T, Patches... &gt; &gt; 

```C++
constexpr bool enable_multipatch_type< MultipatchType< T, Patches... > >;
```




<hr>



### variable is\_multipatch\_type\_v 

```C++
constexpr bool is_multipatch_type_v;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/data_types/multipatch_type.hpp`

