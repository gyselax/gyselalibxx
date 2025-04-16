

# Class IPartialDerivative

**template &lt;class IdxRangeFull, class DerivativeDimension&gt;**



[**ClassList**](annotated.md) **>** [**IPartialDerivative**](classIPartialDerivative.md)



_An abstract class for a partial derivative operator._ [More...](#detailed-description)

* `#include <ipartial_derivative.hpp>`





Inherited by the following classes: [CentralFDMPartialDerivative](classCentralFDMPartialDerivative.md),  [CentralFDMPartialDerivativeWithBValue](classCentralFDMPartialDerivativeWithBValue.md)












## Public Types

| Type | Name |
| ---: | :--- |
| typedef DConstField&lt; IdxRangeFull &gt; | [**DConstFieldType**](#typedef-dconstfieldtype)  <br>_The type of a constant reference to the field to be differentiated._  |
| typedef DField&lt; IdxRangeFull &gt; | [**DFieldType**](#typedef-dfieldtype)  <br>_The type of a reference to the field to be differentiated._  |
| typedef find\_grid\_t&lt; DerivativeDimension, ddc::to\_type\_seq\_t&lt; IdxRangeFull &gt; &gt; | [**GridDerivativeDimension**](#typedef-gridderivativedimension)  <br>_The type of the grid on the dimension on which the partial derivative is calculated._  |
| typedef ddc::remove\_dims\_of\_t&lt; IdxRangeFull, [**GridDerivativeDimension**](classIPartialDerivative.md#typedef-gridderivativedimension) &gt; | [**IdxRangeBatch**](#typedef-idxrangebatch)  <br>_The index range of all dimensions except DerivativeDimension._  |
| typedef IdxRange&lt; [**GridDerivativeDimension**](classIPartialDerivative.md#typedef-gridderivativedimension) &gt; | [**IdxRangeDeriv**](#typedef-idxrangederiv)  <br>_The index range of the dimension on which the partial derivative is calculated._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
| virtual void | [**operator()**](#function-operator) ([**DFieldType**](classIPartialDerivative.md#typedef-dfieldtype) differentiated\_field) const = 0<br>_Compute the partial derivative of a field in a given direction._  |




























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



### typedef GridDerivativeDimension 

_The type of the grid on the dimension on which the partial derivative is calculated._ 
```C++
using IPartialDerivative< IdxRangeFull, DerivativeDimension >::GridDerivativeDimension =  find_grid_t<DerivativeDimension, ddc::to_type_seq_t<IdxRangeFull> >;
```




<hr>



### typedef IdxRangeBatch 

_The index range of all dimensions except DerivativeDimension._ 
```C++
using IPartialDerivative< IdxRangeFull, DerivativeDimension >::IdxRangeBatch =  ddc::remove_dims_of_t<IdxRangeFull, GridDerivativeDimension>;
```




<hr>



### typedef IdxRangeDeriv 

_The index range of the dimension on which the partial derivative is calculated._ 
```C++
using IPartialDerivative< IdxRangeFull, DerivativeDimension >::IdxRangeDeriv =  IdxRange<GridDerivativeDimension>;
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

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/math_tools/ipartial_derivative.hpp`

