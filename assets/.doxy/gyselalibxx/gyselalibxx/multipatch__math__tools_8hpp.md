

# File multipatch\_math\_tools.hpp



[**FileList**](files.md) **>** [**multipatch**](dir_7740c6927b2da0a836b00bedb040a06d.md) **>** [**utils**](dir_573def5310cd01d120c251a7885d602c.md) **>** [**multipatch\_math\_tools.hpp**](multipatch__math__tools_8hpp.md)

[Go to the source code of this file](multipatch__math__tools_8hpp_source.md)



* `#include "l_norm_tools.hpp"`
* `#include "multipatch_field.hpp"`





































## Public Functions

| Type | Name |
| ---: | :--- |
|  double | [**error\_norm\_inf**](#function-error_norm_inf) (ExecSpace exec\_space, detail::MultipatchField&lt; ddc::detail::TypeSeq&lt; Patches... &gt;, TypeSeqInternalTypes &gt; multipatch\_function, detail::MultipatchField&lt; ddc::detail::TypeSeq&lt; Patches... &gt;, TypeSeqInternalTypes &gt; multipatch\_exact\_function) <br>_Compute the infinity norm of the error between 2 Fields or VectorFields over multiple patches._  |
|  double | [**norm\_inf**](#function-norm_inf) (ExecSpace exec\_space, detail::MultipatchField&lt; ddc::detail::TypeSeq&lt; Patches... &gt;, TypeSeqInternalTypes &gt; multipatch\_function) <br>_Compute the infinity norm for a Field or_ [_**VectorField**_](classVectorField.md) _over multiple patches._ |




























## Public Functions Documentation




### function error\_norm\_inf 

_Compute the infinity norm of the error between 2 Fields or VectorFields over multiple patches._ 
```C++
template<class ExecSpace, class... Patches, class TypeSeqInternalTypes>
double error_norm_inf (
    ExecSpace exec_space,
    detail::MultipatchField< ddc::detail::TypeSeq< Patches... >, TypeSeqInternalTypes > multipatch_function,
    detail::MultipatchField< ddc::detail::TypeSeq< Patches... >, TypeSeqInternalTypes > multipatch_exact_function
) 
```





**Parameters:**


* `exec_space` The space on which the function is executed (CPU/GPU). 
* `multipatch_function` The calculated function. 
* `multipatch_exact_function` The exact function with which the calculated function is compared. 



**Returns:**

A double containing the value of the infinity norm. 





        

<hr>



### function norm\_inf 

_Compute the infinity norm for a Field or_ [_**VectorField**_](classVectorField.md) _over multiple patches._
```C++
template<class ExecSpace, class... Patches, class TypeSeqInternalTypes>
double norm_inf (
    ExecSpace exec_space,
    detail::MultipatchField< ddc::detail::TypeSeq< Patches... >, TypeSeqInternalTypes > multipatch_function
) 
```





**Parameters:**


* `exec_space` The space on which the function is executed (CPU/GPU). 
* `multipatch_function` The function whose norm is calculated. 



**Returns:**

A double containing the value of the infinity norm. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/utils/multipatch_math_tools.hpp`

