

# Class IPartialDerivative

**template &lt;class IdxRangeFull, class DerivativeDimension&gt;**



[**ClassList**](annotated.md) **>** [**IPartialDerivative**](classIPartialDerivative.md)



_An abstract class for a partial derivative operator._ [More...](#detailed-description)

* `#include <ipartial_derivative.hpp>`





Inherited by the following classes: [INonAnalyticalPartialDerivative](classINonAnalyticalPartialDerivative.md),  [ConstantPartialDerivative](classConstantPartialDerivative.md),  [INonAnalyticalPartialDerivative](classINonAnalyticalPartialDerivative.md)












## Public Types

| Type | Name |
| ---: | :--- |
| typedef DConstField&lt; IdxRangeFull &gt; | [**DConstFieldType**](#typedef-dconstfieldtype)  <br>_The type of a constant reference to the field to be differentiated._  |
| typedef DField&lt; IdxRangeFull &gt; | [**DFieldType**](#typedef-dfieldtype)  <br>_The type of a reference to the field to be differentiated._  |
| typedef DerivativeDimension | [**DerivDim**](#typedef-derivdim)  <br>_The dimension on which the partial derivative is calculated._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
| virtual void | [**operator()**](#function-operator) ([**DFieldType**](classIPartialDerivative.md#typedef-dfieldtype) differentiated\_field) const = 0<br>_Compute the partial derivative of a field in a given direction._  |
| virtual  | [**~IPartialDerivative**](#function-ipartialderivative) () = default<br> |




























## Detailed Description




**Template parameters:**


* `IdxRangeFull` The index range of the field on which the operator acts (with all dimensions, batched and dimension of interest). 
* `DerivativeDimension` The dimension on which the partial derivative is calculated. 




    
## Public Types Documentation




### typedef DConstFieldType 

_The type of a constant reference to the field to be differentiated._ 
```C++
using IPartialDerivative< IdxRangeFull, DerivativeDimension >::DConstFieldType =  DConstField<IdxRangeFull>;
```




<hr>



### typedef DFieldType 

_The type of a reference to the field to be differentiated._ 
```C++
using IPartialDerivative< IdxRangeFull, DerivativeDimension >::DFieldType =  DField<IdxRangeFull>;
```




<hr>



### typedef DerivDim 

_The dimension on which the partial derivative is calculated._ 
```C++
using IPartialDerivative< IdxRangeFull, DerivativeDimension >::DerivDim =  DerivativeDimension;
```




<hr>
## Public Functions Documentation




### function operator() 

_Compute the partial derivative of a field in a given direction._ 
```C++
virtual void IPartialDerivative::operator() (
    DFieldType differentiated_field
) const = 0
```





**Parameters:**


* `differentiated_field` On output, contains values of the differentiated field. 




        

<hr>



### function ~IPartialDerivative 

```C++
virtual IPartialDerivative::~IPartialDerivative () = default
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/math_tools/ipartial_derivative.hpp`

