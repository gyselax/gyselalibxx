

# Class OrthogonalCoordTransforms

**template &lt;class ArgCoord, class ResultCoord, class JacobianCoord, class... CoordTransform&gt;**



[**ClassList**](annotated.md) **>** [**OrthogonalCoordTransforms**](classOrthogonalCoordTransforms.md)



_A multi-dimensional coordinate transformation which can be decomposed into multiple orthogonal coordinate transformations._ [More...](#detailed-description)

* `#include <orthogonal_coord_transforms.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef ArgCoord | [**CoordArg**](#typedef-coordarg)  <br>_The type of the argument of the function described by this mapping._  |
| typedef JacobianCoord | [**CoordJacobian**](#typedef-coordjacobian)  <br>_The type of the coordinate that can be used to evaluate the Jacobian of this mapping._  |
| typedef ResultCoord | [**CoordResult**](#typedef-coordresult)  <br>_The type of the result of the function described by this mapping._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
|  KOKKOS\_FUNCTION | [**OrthogonalCoordTransforms**](#function-orthogonalcoordtransforms) (CoordTransform const &... transform) <br>_Construct a multi-dimensional coordinate transformation._  |
|  KOKKOS\_INLINE\_FUNCTION auto | [**get\_inverse\_mapping**](#function-get_inverse_mapping) () const<br>_Get the inverse mapping._  |
|  KOKKOS\_FUNCTION double | [**jacobian**](#function-jacobian) ([**CoordJacobian**](classOrthogonalCoordTransforms.md#typedef-coordjacobian) const & coord) const<br>_Compute the Jacobian, the determinant of the Jacobian matrix of the mapping._  |
|  KOKKOS\_FUNCTION double | [**jacobian\_component**](#function-jacobian_component) ([**CoordArg**](classOrthogonalCoordTransforms.md#typedef-coordarg) const & coord) const<br>_Compute the (i,j) coefficient of the Jacobian matrix._  |
|  KOKKOS\_FUNCTION [**DTensor**](classTensor.md)&lt; ddc::to\_type\_seq\_t&lt; [**CoordResult**](classOrthogonalCoordTransforms.md#typedef-coordresult) &gt;, get\_covariant\_dims\_t&lt; ddc::to\_type\_seq\_t&lt; [**CoordArg**](classOrthogonalCoordTransforms.md#typedef-coordarg) &gt; &gt; &gt; | [**jacobian\_matrix**](#function-jacobian_matrix) ([**CoordArg**](classOrthogonalCoordTransforms.md#typedef-coordarg) const & coord) const<br>_Compute full Jacobian matrix._  |
|  KOKKOS\_FUNCTION [**CoordResult**](classOrthogonalCoordTransforms.md#typedef-coordresult) | [**operator()**](#function-operator) ([**CoordArg**](classOrthogonalCoordTransforms.md#typedef-coordarg) const & coord) const<br>_Convert the coordinate to the output coordinate system._  |
|  KOKKOS\_INLINE\_FUNCTION auto | [**operator()**](#function-operator_1) (CoordType const & coord) const<br>_Convert a subset of coordinates to the corresponding subset of the output coordinate system._  |




























## Detailed Description


E.g. \((x_1,x_2) = (y_1 + 3, 4*y_2 + 7)\) Here the coordinates \(y_1\) and \(y_2\) only appear in the description of either \(x_1\) or \(x_2\).


E.g. a cylindrical transformation which can be decomposed into a 2D circular transformation and a linear transformation.




**Template parameters:**


* `ArgCoord` The type of the input coordinates. 
* `ResultCoord` The type of the output coordinates. 
* `JacobianCoord` The type of the coordinates used to calculate the Jacobian. Usually this is the same as ArgCoord. 
* `CoordTransform` The coordinate transformations comprising this coordinate transformation. Note that the order of these is unrelated to the ordering chosen for the coordinates. 




    
## Public Types Documentation




### typedef CoordArg 

_The type of the argument of the function described by this mapping._ 
```C++
using OrthogonalCoordTransforms< ArgCoord, ResultCoord, JacobianCoord, CoordTransform >::CoordArg =  ArgCoord;
```




<hr>



### typedef CoordJacobian 

_The type of the coordinate that can be used to evaluate the Jacobian of this mapping._ 
```C++
using OrthogonalCoordTransforms< ArgCoord, ResultCoord, JacobianCoord, CoordTransform >::CoordJacobian =  JacobianCoord;
```




<hr>



### typedef CoordResult 

_The type of the result of the function described by this mapping._ 
```C++
using OrthogonalCoordTransforms< ArgCoord, ResultCoord, JacobianCoord, CoordTransform >::CoordResult =  ResultCoord;
```




<hr>
## Public Functions Documentation




### function OrthogonalCoordTransforms 

_Construct a multi-dimensional coordinate transformation._ 
```C++
inline explicit KOKKOS_FUNCTION OrthogonalCoordTransforms::OrthogonalCoordTransforms (
    CoordTransform const &... transform
) 
```





**Parameters:**


* `transform` The coordinate transformations comprising this coordinate transformation. 




        

<hr>



### function get\_inverse\_mapping 

_Get the inverse mapping._ 
```C++
inline KOKKOS_INLINE_FUNCTION auto OrthogonalCoordTransforms::get_inverse_mapping () const
```





**Returns:**

The inverse mapping. 





        

<hr>



### function jacobian 

_Compute the Jacobian, the determinant of the Jacobian matrix of the mapping._ 
```C++
inline KOKKOS_FUNCTION double OrthogonalCoordTransforms::jacobian (
    CoordJacobian const & coord
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
inline KOKKOS_FUNCTION double OrthogonalCoordTransforms::jacobian_component (
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
inline KOKKOS_FUNCTION DTensor < ddc::to_type_seq_t< CoordResult >, get_covariant_dims_t< ddc::to_type_seq_t< CoordArg > > > OrthogonalCoordTransforms::jacobian_matrix (
    CoordArg const & coord
) const
```



For some computations, we need the complete Jacobian matrix or just the coefficients. The coefficients can be given independently with the function jacobian\_component.




**Parameters:**


* `coord` The coordinate where we evaluate the Jacobian matrix. 



**Returns:**

The Jacobian matrix. 





        

<hr>



### function operator() 

_Convert the coordinate to the output coordinate system._ 
```C++
inline KOKKOS_FUNCTION CoordResult OrthogonalCoordTransforms::operator() (
    CoordArg const & coord
) const
```





**Parameters:**


* `coord` The coordinate to be converted expressed on the input coordinate system.



**Returns:**

The coordinate expressed on the output coordinate system. 





        

<hr>



### function operator() 

_Convert a subset of coordinates to the corresponding subset of the output coordinate system._ 
```C++
template<class CoordType, std::enable_if_t<!std::is_same_v< CoordType, CoordArg >, bool >>
inline KOKKOS_INLINE_FUNCTION auto OrthogonalCoordTransforms::operator() (
    CoordType const & coord
) const
```



Neglected elements of the coordinate system must be orthogonal to the provided coordinates.




**Parameters:**


* `coord` The coordinate to be converted expressed on the input coordinate system.



**Returns:**

The coordinate expressed on the output coordinate system. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/coord_transformations/orthogonal_coord_transforms.hpp`

