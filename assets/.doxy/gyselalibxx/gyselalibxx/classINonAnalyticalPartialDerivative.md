

# Class INonAnalyticalPartialDerivative

**template &lt;class IdxRangeFull, class DerivativeDimension&gt;**



[**ClassList**](annotated.md) **>** [**INonAnalyticalPartialDerivative**](classINonAnalyticalPartialDerivative.md)



_An abstract class for a non-analytical partial derivative operator. This class inherits from_ [_**IPartialDerivative**_](classIPartialDerivative.md) _but adds IdxRange information which is useful when iterating over the dimension along which the derivative is defined. This information may not be relevant for an analytical_[_**IPartialDerivative**_](classIPartialDerivative.md) _as in this case no IdxRange iteration is required._[More...](#detailed-description)

* `#include <ipartial_derivative.hpp>`



Inherits the following classes: [IPartialDerivative](classIPartialDerivative.md)


Inherited by the following classes: [CentralFDMPartialDerivative](classCentralFDMPartialDerivative.md),  [CentralFDMPartialDerivativeWithBValue](classCentralFDMPartialDerivativeWithBValue.md)












## Public Types

| Type | Name |
| ---: | :--- |
| typedef find\_grid\_t&lt; DerivativeDimension, ddc::to\_type\_seq\_t&lt; IdxRangeFull &gt; &gt; | [**GridDerivativeDimension**](#typedef-gridderivativedimension)  <br>_The type of the grid on the dimension on which the partial derivative is calculated._  |
| typedef ddc::remove\_dims\_of\_t&lt; IdxRangeFull, [**GridDerivativeDimension**](classINonAnalyticalPartialDerivative.md#typedef-gridderivativedimension) &gt; | [**IdxRangeBatch**](#typedef-idxrangebatch)  <br>_The index range of all dimensions except DerivativeDimension._  |
| typedef IdxRange&lt; [**GridDerivativeDimension**](classINonAnalyticalPartialDerivative.md#typedef-gridderivativedimension) &gt; | [**IdxRangeDeriv**](#typedef-idxrangederiv)  <br>_The index range of the dimension on which the partial derivative is calculated._  |


## Public Types inherited from IPartialDerivative

See [IPartialDerivative](classIPartialDerivative.md)

| Type | Name |
| ---: | :--- |
| typedef DConstField&lt; IdxRangeFull &gt; | [**DConstFieldType**](classIPartialDerivative.md#typedef-dconstfieldtype)  <br>_The type of a constant reference to the field to be differentiated._  |
| typedef DField&lt; IdxRangeFull &gt; | [**DFieldType**](classIPartialDerivative.md#typedef-dfieldtype)  <br>_The type of a reference to the field to be differentiated._  |
| typedef DerivativeDimension | [**DerivDim**](classIPartialDerivative.md#typedef-derivdim)  <br>_The dimension on which the partial derivative is calculated._  |








































## Public Functions inherited from IPartialDerivative

See [IPartialDerivative](classIPartialDerivative.md)

| Type | Name |
| ---: | :--- |
| virtual void | [**operator()**](classIPartialDerivative.md#function-operator) ([**DFieldType**](classIPartialDerivative.md#typedef-dfieldtype) differentiated\_field) const = 0<br>_Compute the partial derivative of a field in a given direction._  |
| virtual  | [**~IPartialDerivative**](classIPartialDerivative.md#function-ipartialderivative) () = default<br> |






















































## Detailed Description




**Template parameters:**


* `IdxRangeFull` The index range of the field on which the operator acts (with all dimensions, batched and dimension of interest). 
* `DerivativeDimension` The dimension on which the partial derivative is calculated. 




    
## Public Types Documentation




### typedef GridDerivativeDimension 

_The type of the grid on the dimension on which the partial derivative is calculated._ 
```C++
using INonAnalyticalPartialDerivative< IdxRangeFull, DerivativeDimension >::GridDerivativeDimension =  find_grid_t<DerivativeDimension, ddc::to_type_seq_t<IdxRangeFull>>;
```




<hr>



### typedef IdxRangeBatch 

_The index range of all dimensions except DerivativeDimension._ 
```C++
using INonAnalyticalPartialDerivative< IdxRangeFull, DerivativeDimension >::IdxRangeBatch =  ddc::remove_dims_of_t<IdxRangeFull, GridDerivativeDimension>;
```




<hr>



### typedef IdxRangeDeriv 

_The index range of the dimension on which the partial derivative is calculated._ 
```C++
using INonAnalyticalPartialDerivative< IdxRangeFull, DerivativeDimension >::IdxRangeDeriv =  IdxRange<GridDerivativeDimension>;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/math_tools/ipartial_derivative.hpp`

