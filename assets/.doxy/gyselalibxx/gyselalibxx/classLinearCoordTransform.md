

# Class LinearCoordTransform

**template &lt;class InputDim, class OutputDim, class CoordJacob&gt;**



[**ClassList**](annotated.md) **>** [**LinearCoordTransform**](classLinearCoordTransform.md)



_A class describing a linear coordinate transformation._ [More...](#detailed-description)

* `#include <linear_coord_transform.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef Coord&lt; InputDim &gt; | [**CoordArg**](#typedef-coordarg)  <br>_The type of the argument of the function described by this mapping._  |
| typedef CoordJacob | [**CoordJacobian**](#typedef-coordjacobian)  <br>_The type of the coordinate that can be used to evaluate the Jacobian of this mapping._  |
| typedef Coord&lt; OutputDim &gt; | [**CoordResult**](#typedef-coordresult)  <br>_The type of the result of the function described by this mapping._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
|  KOKKOS\_FUNCTION | [**LinearCoordTransform**](#function-linearcoordtransform-12) (Coord&lt; InputDim &gt; reference\_point\_on\_input\_dim, Coord&lt; OutputDim &gt; reference\_point\_on\_output\_dim, double scaling\_factor) <br>_A constructor for the linear coordinate transformation._  |
|  KOKKOS\_DEFAULTED\_FUNCTION | [**LinearCoordTransform**](#function-linearcoordtransform-22) ([**LinearCoordTransform**](classLinearCoordTransform.md) const & other) = default<br> |
|  KOKKOS\_INLINE\_FUNCTION [**LinearCoordTransform**](classLinearCoordTransform.md)&lt; OutputDim, InputDim &gt; | [**get\_inverse\_mapping**](#function-get_inverse_mapping) () const<br>_Get the inverse mapping._  |
|  KOKKOS\_FUNCTION double | [**inv\_jacobian\_component**](#function-inv_jacobian_component) ([**CoordArg**](classLinearCoordTransform.md#typedef-coordarg) const & coord) const<br>_Compute the (i,j) coefficient of the inverse Jacobian matrix._  |
|  KOKKOS\_FUNCTION [**DTensor**](classTensor.md)&lt; VectorIndexSet&lt; InputDim &gt;, VectorIndexSet&lt; OutputDim &gt; &gt; | [**inv\_jacobian\_matrix**](#function-inv_jacobian_matrix) ([**CoordArg**](classLinearCoordTransform.md#typedef-coordarg) const & coord) const<br>_Compute full inverse Jacobian matrix._  |
|  KOKKOS\_FUNCTION double | [**jacobian**](#function-jacobian) ([**CoordArg**](classLinearCoordTransform.md#typedef-coordarg) const & coord) const<br>_Compute the Jacobian, the determinant of the Jacobian matrix of the mapping._  |
|  KOKKOS\_FUNCTION double | [**jacobian\_component**](#function-jacobian_component) ([**CoordArg**](classLinearCoordTransform.md#typedef-coordarg) const & coord) const<br>_Compute the (i,j) coefficient of the Jacobian matrix._  |
|  KOKKOS\_FUNCTION [**DTensor**](classTensor.md)&lt; VectorIndexSet&lt; OutputDim &gt;, VectorIndexSet&lt; InputDim &gt; &gt; | [**jacobian\_matrix**](#function-jacobian_matrix) ([**CoordArg**](classLinearCoordTransform.md#typedef-coordarg) const & coord) const<br>_Compute full Jacobian matrix._  |
|  KOKKOS\_FUNCTION [**CoordResult**](classLinearCoordTransform.md#typedef-coordresult) | [**operator()**](#function-operator) ([**CoordArg**](classLinearCoordTransform.md#typedef-coordarg) const & coord) const<br>_Convert the coordinate on the input dimension to a coordinate on the output dimension._  |
|  KOKKOS\_DEFAULTED\_FUNCTION | [**~LinearCoordTransform**](#function-linearcoordtransform) () = default<br> |




























## Detailed Description


A class describing a linear coordinate transformation of the form: \(x_1 = \alpha x_2 + \beta\) where \(x_1\) and \(x_2\) are coordinates in the input and output dimensions, and \(\alpha\) and \(\beta\) are coefficients. 


    
## Public Types Documentation




### typedef CoordArg 

_The type of the argument of the function described by this mapping._ 
```C++
using LinearCoordTransform< InputDim, OutputDim, CoordJacob >::CoordArg =  Coord<InputDim>;
```




<hr>



### typedef CoordJacobian 

_The type of the coordinate that can be used to evaluate the Jacobian of this mapping._ 
```C++
using LinearCoordTransform< InputDim, OutputDim, CoordJacob >::CoordJacobian =  CoordJacob;
```




<hr>



### typedef CoordResult 

_The type of the result of the function described by this mapping._ 
```C++
using LinearCoordTransform< InputDim, OutputDim, CoordJacob >::CoordResult =  Coord<OutputDim>;
```




<hr>
## Public Functions Documentation




### function LinearCoordTransform [1/2]

_A constructor for the linear coordinate transformation._ 
```C++
inline explicit KOKKOS_FUNCTION LinearCoordTransform::LinearCoordTransform (
    Coord< InputDim > reference_point_on_input_dim,
    Coord< OutputDim > reference_point_on_output_dim,
    double scaling_factor
) 
```



The constructor infers the linear coordinate transformation from a reference point provided in both coordinate systems and a scaling factor.


The transformation is then defined as: \(x_{out} = x_{out}^* + s * (x_{in} - x_{in}^*)\) where \(x_{in}\) and \(x_{out}\) denote the point expressed in the input and output coordinate system, \(\cdot^*\) denotes the reference point and \(s\) denotes the scaling factor.




**Parameters:**


* `reference_point_on_input_dim` The reference point expressed in the input coordinate system. 
* `reference_point_on_output_dim` The reference point expressed in the output coordinate system. 
* `scaling_factor` The scaling factor describing how distances in the output coordinate system scale compared to distances in the input coordinate system. 




        

<hr>



### function LinearCoordTransform [2/2]

```C++
KOKKOS_DEFAULTED_FUNCTION LinearCoordTransform::LinearCoordTransform (
    LinearCoordTransform const & other
) = default
```



Copy constructor for [**LinearCoordTransform**](classLinearCoordTransform.md). 

**Parameters:**


* `other` The operator to be copied. 




        

<hr>



### function get\_inverse\_mapping 

_Get the inverse mapping._ 
```C++
inline KOKKOS_INLINE_FUNCTION LinearCoordTransform < OutputDim, InputDim > LinearCoordTransform::get_inverse_mapping () const
```





**Returns:**

The inverse mapping. 





        

<hr>



### function inv\_jacobian\_component 

_Compute the (i,j) coefficient of the inverse Jacobian matrix._ 
```C++
template<class IndexTag1, class IndexTag2>
inline KOKKOS_FUNCTION double LinearCoordTransform::inv_jacobian_component (
    CoordArg const & coord
) const
```





**Parameters:**


* `coord` The coordinate where we evaluate the inverse Jacobian matrix.



**Returns:**

A double with the value of the (i,j) coefficient of the inverse Jacobian matrix. 





        

<hr>



### function inv\_jacobian\_matrix 

_Compute full inverse Jacobian matrix._ 
```C++
inline KOKKOS_FUNCTION DTensor < VectorIndexSet< InputDim >, VectorIndexSet< OutputDim > > LinearCoordTransform::inv_jacobian_matrix (
    CoordArg const & coord
) const
```



For some computations, we need the complete inverse Jacobian matrix or just the coefficients. The coefficient can be given as a scalar with the function inv\_jacobian\_component.




**Parameters:**


* `coord` The coordinate where we evaluate the Jacobian matrix. 



**Returns:**

The inverse Jacobian matrix. 





        

<hr>



### function jacobian 

_Compute the Jacobian, the determinant of the Jacobian matrix of the mapping._ 
```C++
inline KOKKOS_FUNCTION double LinearCoordTransform::jacobian (
    CoordArg const & coord
) const
```





**Parameters:**


* `coord` The coordinate where we evaluate the Jacobian.



**Returns:**

A double with the value of the determinant of the Jacobian matrix. 





        

<hr>



### function jacobian\_component 

_Compute the (i,j) coefficient of the Jacobian matrix._ 
```C++
template<class IndexTag1, class IndexTag2>
inline KOKKOS_FUNCTION double LinearCoordTransform::jacobian_component (
    CoordArg const & coord
) const
```





**Parameters:**


* `coord` The coordinate where we evaluate the Jacobian matrix.



**Returns:**

A double with the value of the (i,j) coefficient of the Jacobian matrix. 





        

<hr>



### function jacobian\_matrix 

_Compute full Jacobian matrix._ 
```C++
inline KOKKOS_FUNCTION DTensor < VectorIndexSet< OutputDim >, VectorIndexSet< InputDim > > LinearCoordTransform::jacobian_matrix (
    CoordArg const & coord
) const
```



For some computations, we need the complete Jacobian matrix or just the coefficients. The coefficient can be given as a scalar with the function jacobian\_component.




**Parameters:**


* `coord` The coordinate where we evaluate the Jacobian matrix. 



**Returns:**

The Jacobian matrix. 





        

<hr>



### function operator() 

_Convert the coordinate on the input dimension to a coordinate on the output dimension._ 
```C++
inline KOKKOS_FUNCTION CoordResult LinearCoordTransform::operator() (
    CoordArg const & coord
) const
```





**Parameters:**


* `coord` The coordinate to be converted expressed on the input coordinate system.



**Returns:**

The coordinate expressed on the output coordinate system. 





        

<hr>



### function ~LinearCoordTransform 

```C++
KOKKOS_DEFAULTED_FUNCTION LinearCoordTransform::~LinearCoordTransform () = default
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/coord_transformations/linear_coord_transform.hpp`

