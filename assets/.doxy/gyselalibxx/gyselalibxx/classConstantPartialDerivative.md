

# Class ConstantPartialDerivative

**template &lt;class IdxRangeFull, class DerivativeDimension&gt;**



[**ClassList**](annotated.md) **>** [**ConstantPartialDerivative**](classConstantPartialDerivative.md)



_A class to get the derivative of a constant function. When the derivative of a function is known to be 0 but the dimension is still needed this class can be used to avoid unnecessary calculations._ [More...](#detailed-description)

* `#include <constant_partial_derivatives.hpp>`



Inherits the following classes: [IPartialDerivative](classIPartialDerivative.md)
















## Public Types inherited from IPartialDerivative

See [IPartialDerivative](classIPartialDerivative.md)

| Type | Name |
| ---: | :--- |
| typedef DConstField&lt; IdxRangeFull &gt; | [**DConstFieldType**](classIPartialDerivative.md#typedef-dconstfieldtype)  <br>_The type of a constant reference to the field to be differentiated._  |
| typedef DField&lt; IdxRangeFull &gt; | [**DFieldType**](classIPartialDerivative.md#typedef-dfieldtype)  <br>_The type of a reference to the field to be differentiated._  |
| typedef find\_grid\_t&lt; DerivativeDimension, ddc::to\_type\_seq\_t&lt; IdxRangeFull &gt; &gt; | [**GridDerivativeDimension**](classIPartialDerivative.md#typedef-gridderivativedimension)  <br>_The type of the grid on the dimension on which the partial derivative is calculated._  |
| typedef ddc::remove\_dims\_of\_t&lt; IdxRangeFull, [**GridDerivativeDimension**](classIPartialDerivative.md#typedef-gridderivativedimension) &gt; | [**IdxRangeBatch**](classIPartialDerivative.md#typedef-idxrangebatch)  <br>_The index range of all dimensions except DerivativeDimension._  |
| typedef IdxRange&lt; [**GridDerivativeDimension**](classIPartialDerivative.md#typedef-gridderivativedimension) &gt; | [**IdxRangeDeriv**](classIPartialDerivative.md#typedef-idxrangederiv)  <br>_The index range of the dimension on which the partial derivative is calculated._  |






































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**ConstantPartialDerivative**](#function-constantpartialderivative) (double deriv\_value) <br>_Create an instance of_ [_**ConstantPartialDerivative**_](classConstantPartialDerivative.md) _._ |
| virtual void | [**operator()**](#function-operator) (DFieldType differentiated\_field) const<br>_Set the partial derivative of a field to 0._  |


## Public Functions inherited from IPartialDerivative

See [IPartialDerivative](classIPartialDerivative.md)

| Type | Name |
| ---: | :--- |
| virtual void | [**operator()**](classIPartialDerivative.md#function-operator) ([**DFieldType**](classIPartialDerivative.md#typedef-dfieldtype) differentiated\_field) const = 0<br>_Compute the partial derivative of a field in a given direction._  |
| virtual  | [**~IPartialDerivative**](classIPartialDerivative.md#function-ipartialderivative) () = default<br> |






















































## Detailed Description




**Template parameters:**


* `IdxRangeFull` The index range of the field on which the operator acts (with all dimensions, batched and dimension of interest, used for inheritance). 
* `DerivativeDimension` The dimension on which the partial derivative is calculated. (used for inheritance). 




    
## Public Functions Documentation




### function ConstantPartialDerivative 

_Create an instance of_ [_**ConstantPartialDerivative**_](classConstantPartialDerivative.md) _._
```C++
inline explicit ConstantPartialDerivative::ConstantPartialDerivative (
    double deriv_value
) 
```





**Parameters:**


* `deriv_value` The value that should be returned as the constant value of the derivative. 




        

<hr>



### function operator() 

_Set the partial derivative of a field to 0._ 
```C++
inline virtual void ConstantPartialDerivative::operator() (
    DFieldType differentiated_field
) const
```





**Parameters:**


* `differentiated_field` On output, contains values of the differentiated field. 




        
Implements [*IPartialDerivative::operator()*](classIPartialDerivative.md#function-operator)


<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/math_tools/constant_partial_derivatives.hpp`

