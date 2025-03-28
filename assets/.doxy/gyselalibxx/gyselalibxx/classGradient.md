

# Class Gradient

**template &lt;class MetricTensorType&gt;**



[**ClassList**](annotated.md) **>** [**Gradient**](classGradient.md)



_A class which implements a gradient operator._ [More...](#detailed-description)

* `#include <gradient.hpp>`





































## Public Functions

| Type | Name |
| ---: | :--- |
|  KOKKOS\_FUNCTION | [**Gradient**](#function-gradient) (MetricTensorType const & metric\_tensor) <br>_Construct an instance of the class_ [_**Gradient**_](classGradient.md) _._ |
|  KOKKOS\_INLINE\_FUNCTION DVectorCov | [**operator()**](#function-operator) (DVectorCov const & partial\_derivatives) const<br>_Compute the gradient of a scalar field at a given coordinate, from the partial derivatives of the field. The gradient is expressed on the covariant basis. The components of the gradient in the covariant basis are simply equal to the value of the partial derivatives of the scalar field._  |
|  KOKKOS\_INLINE\_FUNCTION void | [**operator()**](#function-operator_1) ([**DVectorFieldCovType**](classVectorField.md)&lt; IdxRange &gt; gradient, [**DVectorConstFieldCovType**](classVectorField.md)&lt; IdxRange &gt; const partial\_derivatives) const<br>_Compute the gradient of a scalar field at a given coordinate, using partial derivatives of the field. The gradient is expressed on the covariant basis. The components of the gradient in the covariant basis are simply equal to the value of the partial derivatives of the scalar field._  |
|  KOKKOS\_INLINE\_FUNCTION DVectorType | [**operator()**](#function-operator_2) (DVectorCov const & partial\_derivatives, CoordArg const & coord) const<br>_Compute the gradient of a scalar field at a given coordinate, using partial derivatives of the field. The gradient is expressed on the contravariant basis._  __ |
|  KOKKOS\_FUNCTION void | [**operator()**](#function-operator_3) ([**DVectorFieldType**](classVectorField.md)&lt; IdxRange &gt; const gradient, [**DVectorConstFieldCovType**](classVectorField.md)&lt; IdxRange &gt; const partial\_derivatives) const<br>_Compute the gradient of a scalar field using partial derivatives of the field. The gradient is expressed on the contravariant basis._  __ |




























## Detailed Description




**Template parameters:**


* `MetricTensorType` A type representing a metric tensor. 




    
## Public Functions Documentation




### function Gradient 

_Construct an instance of the class_ [_**Gradient**_](classGradient.md) _._
```C++
inline explicit KOKKOS_FUNCTION Gradient::Gradient (
    MetricTensorType const & metric_tensor
) 
```





**Parameters:**


* `metric_tensor` A [**MetricTensorEvaluator**](classMetricTensorEvaluator.md). 




        

<hr>



### function operator() 

_Compute the gradient of a scalar field at a given coordinate, from the partial derivatives of the field. The gradient is expressed on the covariant basis. The components of the gradient in the covariant basis are simply equal to the value of the partial derivatives of the scalar field._ 
```C++
inline KOKKOS_INLINE_FUNCTION DVectorCov Gradient::operator() (
    DVectorCov const & partial_derivatives
) const
```





**Parameters:**


* `partial_derivatives` A vector containing the partial derivatives of the scalar field expressed at a given coordinate.



**Returns:**

The components of the gradient at a given coordinate, expressed on the covariant basis. 





        

<hr>



### function operator() 

_Compute the gradient of a scalar field at a given coordinate, using partial derivatives of the field. The gradient is expressed on the covariant basis. The components of the gradient in the covariant basis are simply equal to the value of the partial derivatives of the scalar field._ 
```C++
template<class IdxRange>
inline KOKKOS_INLINE_FUNCTION void Gradient::operator() (
    DVectorFieldCovType < IdxRange > gradient,
    DVectorConstFieldCovType < IdxRange > const partial_derivatives
) const
```





**Parameters:**


* `gradient` A vector field that contains on output the value of the gradient components expressed on the contravariant 
* `partial_derivatives` A vector field containing the partial derivatives of the scalar field. basis. 




        

<hr>



### function operator() 

_Compute the gradient of a scalar field at a given coordinate, using partial derivatives of the field. The gradient is expressed on the contravariant basis._  __
```C++
inline KOKKOS_INLINE_FUNCTION DVectorType Gradient::operator() (
    DVectorCov const & partial_derivatives,
    CoordArg const & coord
) const
```





**Parameters:**


* `partial_derivatives` A vector that contains the partial derivatives of the scalar field expressed at a given coordinate. 
* `coord` The coordinate at which the gradient should be evaluated.



**Returns:**

The components of the gradient at a given coordinate, expressed on the contravariant basis. 





        

<hr>



### function operator() 

_Compute the gradient of a scalar field using partial derivatives of the field. The gradient is expressed on the contravariant basis._  __
```C++
template<class IdxRange>
inline KOKKOS_FUNCTION void Gradient::operator() (
    DVectorFieldType < IdxRange > const gradient,
    DVectorConstFieldCovType < IdxRange > const partial_derivatives
) const
```





**Parameters:**


* `gradient` A vector field that contains on output the value of the gradient components expressed on the contravariant 
* `partial_derivatives` A vector field that contains the partial derivatives of the scalar field. basis. 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/math_tools/gradient.hpp`

