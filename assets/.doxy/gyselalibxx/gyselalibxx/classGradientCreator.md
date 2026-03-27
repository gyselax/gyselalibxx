

# Class GradientCreator

**template &lt;typename IdxRangeFull, typename... DerivativeDims&gt;**



[**ClassList**](annotated.md) **>** [**GradientCreator**](classGradientCreator.md)



_Operator to calculate the gradient of a function._ [More...](#detailed-description)

* `#include <gradient_creator.hpp>`





































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**GradientCreator**](#function-gradientcreator) ([**IPartialDerivativeCreator**](classIPartialDerivativeCreator.md)&lt; IdxRangeFull, DerivativeDims &gt; const &... partial\_derivative\_operator) <br>_Create a_ [_**GradientCreator**_](classGradientCreator.md) _._ |
|  void | [**operator()**](#function-operator) ([**DVectorField**](classVectorField.md)&lt; IdxRangeFull, get\_covariant\_dims\_t&lt; VectorIndexSet&lt; DerivativeDims... &gt; &gt; &gt; grad\_func\_cov, DConstField&lt; IdxRangeFull &gt; func) const<br>_Fill a_ [_**VectorField**_](classVectorField.md) _with the values of the gradient of a function._ |




























## Detailed Description




**Template parameters:**


* `IdxRangeFull` The index range on which the function being differentiated should be defined. 
* `DerivativeDims` The dimensions along which the function should be differentiated to construct the gradient. 




    
## Public Functions Documentation




### function GradientCreator 

_Create a_ [_**GradientCreator**_](classGradientCreator.md) _._
```C++
inline explicit GradientCreator::GradientCreator (
    IPartialDerivativeCreator < IdxRangeFull, DerivativeDims > const &... partial_derivative_operator
) 
```





**Parameters:**


* `partial_derivative_operator` The operators used to calculate the derivative along each of the dimensions. 




        

<hr>



### function operator() 

_Fill a_ [_**VectorField**_](classVectorField.md) _with the values of the gradient of a function._
```C++
inline void GradientCreator::operator() (
    DVectorField < IdxRangeFull, get_covariant_dims_t< VectorIndexSet< DerivativeDims... > > > grad_func_cov,
    DConstField< IdxRangeFull > func
) const
```





**Parameters:**


* `grad_func_cov` The gradient of the function expressed as a covariant vector. 
* `func` The function whose gradient is calculated. 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/math_tools/gradient_creator.hpp`

