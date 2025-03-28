

# File derivative\_field\_common.hpp



[**FileList**](files.md) **>** [**data\_types**](dir_eaa769653453aaefd8cc10e98e9bb3eb.md) **>** [**derivative\_field\_common.hpp**](derivative__field__common_8hpp.md)

[Go to the source code of this file](derivative__field__common_8hpp_source.md)



* `#include <array>`
* `#include <ddc/ddc.hpp>`
* `#include <ddc/kernels/splines.hpp>`
* `#include "ddc_aliases.hpp"`
* `#include "deriv_details.hpp"`
* `#include "idx_range_slice.hpp"`















## Classes

| Type | Name |
| ---: | :--- |
| class | [**DerivFieldCommon&lt; FieldType, IdxRange&lt; DDims... &gt; &gt;**](classDerivFieldCommon_3_01FieldType_00_01IdxRange_3_01DDims_8_8_8_01_4_01_4.md) &lt;class FieldType, DDims&gt;<br>_An abstract class which holds a chunk of memory describing a field and its derivatives. This is the superclass for_ [_**DerivFieldMem**_](classDerivFieldMem.md) _and_[_**DerivField**_](classDerivField.md) _._ |






## Public Attributes

| Type | Name |
| ---: | :--- |
|  constexpr bool | [**enable\_borrowed\_deriv\_field**](#variable-enable_borrowed_deriv_field)   = `false`<br> |
|  constexpr bool | [**enable\_deriv\_field**](#variable-enable_deriv_field)   = `false`<br> |
|  constexpr bool | [**is\_borrowed\_deriv\_field\_v**](#variable-is_borrowed_deriv_field_v)   = `/* multi line expression */`<br> |
|  constexpr bool | [**is\_deriv\_field\_v**](#variable-is_deriv_field_v)   = `enable\_deriv\_field&lt;std::remove\_const\_t&lt;std::remove\_reference\_t&lt;[**T**](structT.md)&gt;&gt;&gt;`<br> |












































## Public Attributes Documentation




### variable enable\_borrowed\_deriv\_field 

```C++
constexpr bool enable_borrowed_deriv_field;
```




<hr>



### variable enable\_deriv\_field 

```C++
constexpr bool enable_deriv_field;
```




<hr>



### variable is\_borrowed\_deriv\_field\_v 

```C++
constexpr bool is_borrowed_deriv_field_v;
```




<hr>



### variable is\_deriv\_field\_v 

```C++
constexpr bool is_deriv_field_v;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/data_types/derivative_field_common.hpp`

