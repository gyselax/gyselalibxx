

# File multipatch\_math\_tools.hpp



[**FileList**](files.md) **>** [**multipatch**](dir_7740c6927b2da0a836b00bedb040a06d.md) **>** [**utils**](dir_573def5310cd01d120c251a7885d602c.md) **>** [**multipatch\_math\_tools.hpp**](multipatch__math__tools_8hpp.md)

[Go to the source code of this file](multipatch__math__tools_8hpp_source.md)



* `#include "l_norm_tools.hpp"`
* `#include "multipatch_field.hpp"`





































## Public Functions

| Type | Name |
| ---: | :--- |
|  double | [**error\_norm\_inf**](#function-error_norm_inf) (ExecSpace exec\_space, [**MultipatchField**](classMultipatchField.md)&lt; [**T**](structT.md), Patches... &gt; multipatch\_function, [**MultipatchField**](classMultipatchField.md)&lt; [**T**](structT.md), Patches... &gt; multipatch\_exact\_function) <br>_Compute the infinity norm of the error between 2 Fields or VectorFields over multiple patches._  |
|  double | [**norm\_inf**](#function-norm_inf) (ExecSpace exec\_space, [**MultipatchField**](classMultipatchField.md)&lt; [**T**](structT.md), Patches... &gt; multipatch\_function) <br>_Compute the infinity norm for a Field or_ [_**VectorField**_](classVectorField.md) _over multiple patches._ |




























## Public Functions Documentation




### function error\_norm\_inf 

_Compute the infinity norm of the error between 2 Fields or VectorFields over multiple patches._ 
```C++
template<class ExecSpace, template< typename P > typename T, class... Patches>
double error_norm_inf (
    ExecSpace exec_space,
    MultipatchField < T , Patches... > multipatch_function,
    MultipatchField < T , Patches... > multipatch_exact_function
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
template<class ExecSpace, template< typename P > typename T, class... Patches>
double norm_inf (
    ExecSpace exec_space,
    MultipatchField < T , Patches... > multipatch_function
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

