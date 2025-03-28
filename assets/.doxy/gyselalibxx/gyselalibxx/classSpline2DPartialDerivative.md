

# Class Spline2DPartialDerivative

**template &lt;class [**SplineBuilder2DCache**](classSplineBuilder2DCache.md), class SplineEvaluator2D, class DerivativeDimension&gt;**



[**ClassList**](annotated.md) **>** [**Spline2DPartialDerivative**](classSpline2DPartialDerivative.md)



_A class which implements a partial derivative operator using a 2d spline interpolation._ [More...](#detailed-description)

* `#include <spline_2d_partial_derivative.hpp>`



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
|   | [**Spline2DPartialDerivative**](#function-spline2dpartialderivative) ([**SplineBuilder2DCache**](classSplineBuilder2DCache.md) & builder\_cache, SplineEvaluator2D const & evaluator, DConstFieldType const field) <br>_Construct an instance of the class_ [_**Spline2DPartialDerivative**_](classSpline2DPartialDerivative.md) _._ |
|  void | [**operator()**](#function-operator) (DFieldType differentiated\_field) const<br>_Compute the partial derivative of a field in the direction where the field is represented using 2d splines._  |


## Public Functions inherited from IPartialDerivative

See [IPartialDerivative](classIPartialDerivative.md)

| Type | Name |
| ---: | :--- |
| virtual void | [**operator()**](classIPartialDerivative.md#function-operator) ([**DFieldType**](classIPartialDerivative.md#typedef-dfieldtype) differentiated\_field) const = 0<br>_Compute the partial derivative of a field in a given direction._  |






















































## Detailed Description




**Template parameters:**


* `SplineBuilder2D` A 2D spline builder. 
* `SplineEvaluator2D` A 2D spline evaluator. 




    
## Public Functions Documentation




### function Spline2DPartialDerivative 

_Construct an instance of the class_ [_**Spline2DPartialDerivative**_](classSpline2DPartialDerivative.md) _._
```C++
inline explicit Spline2DPartialDerivative::Spline2DPartialDerivative (
    SplineBuilder2DCache & builder_cache,
    SplineEvaluator2D const & evaluator,
    DConstFieldType const field
) 
```





**Parameters:**


* `builder_cache` A 2D spline builder cache. 
* `evaluator` A 2D spline evaluator. 
* `field` The field to be differentiated. 




        

<hr>



### function operator() 

_Compute the partial derivative of a field in the direction where the field is represented using 2d splines._ 
```C++
inline void Spline2DPartialDerivative::operator() (
    DFieldType differentiated_field
) const
```





**Parameters:**


* `differentiated_field` Contains on output the value of the differentiated field. 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/math_tools/spline_2d_partial_derivative.hpp`

