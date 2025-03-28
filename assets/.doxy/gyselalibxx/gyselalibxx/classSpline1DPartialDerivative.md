

# Class Spline1DPartialDerivative

**template &lt;class Spline1DBuilder, class Spline1DEvaluator&gt;**



[**ClassList**](annotated.md) **>** [**Spline1DPartialDerivative**](classSpline1DPartialDerivative.md)



_A class which implements a partial derivative operator using a 1d spline interpolation._ [More...](#detailed-description)

* `#include <spline_1d_partial_derivative.hpp>`



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
|   | [**Spline1DPartialDerivative**](#function-spline1dpartialderivative) (Spline1DBuilder const & builder, Spline1DEvaluator const & evaluator, DConstFieldType const field) <br>_Construct an instance of the class_ [_**Spline1DPartialDerivative**_](classSpline1DPartialDerivative.md) _._ |
|  void | [**operator()**](#function-operator) (DFieldType differentiated\_field) const<br>_Compute the partial derivative of a field in the direction where the field is represented using 1d splines._  |


## Public Functions inherited from IPartialDerivative

See [IPartialDerivative](classIPartialDerivative.md)

| Type | Name |
| ---: | :--- |
| virtual void | [**operator()**](classIPartialDerivative.md#function-operator) ([**DFieldType**](classIPartialDerivative.md#typedef-dfieldtype) differentiated\_field) const = 0<br>_Compute the partial derivative of a field in a given direction._  |






















































## Detailed Description




**Template parameters:**


* `Spline1DBuilder` A 1D spline builder. 
* `Spline1DEvaluator` A 1D spline evaluator. 




    
## Public Functions Documentation




### function Spline1DPartialDerivative 

_Construct an instance of the class_ [_**Spline1DPartialDerivative**_](classSpline1DPartialDerivative.md) _._
```C++
inline explicit Spline1DPartialDerivative::Spline1DPartialDerivative (
    Spline1DBuilder const & builder,
    Spline1DEvaluator const & evaluator,
    DConstFieldType const field
) 
```





**Parameters:**


* `builder` A 1D spline builder. 
* `evaluator` A 1D spline evaluator. 
* `field` The field to be differentiated. 




        

<hr>



### function operator() 

_Compute the partial derivative of a field in the direction where the field is represented using 1d splines._ 
```C++
inline void Spline1DPartialDerivative::operator() (
    DFieldType differentiated_field
) const
```





**Parameters:**


* `differentiated_field` Contains on output the value of the differentiated field. 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/math_tools/spline_1d_partial_derivative.hpp`

