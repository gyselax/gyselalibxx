

# Class CartesianToBarycentric

**template &lt;class [**X**](structX.md), class [**Y**](structY.md), class Corner1Tag, class Corner2Tag, class Corner3Tag&gt;**



[**ClassList**](annotated.md) **>** [**CartesianToBarycentric**](classCartesianToBarycentric.md)



_A class to convert Cartesian coordinates to barycentric coordinates on a triangle._ [More...](#detailed-description)

* `#include <cartesian_to_barycentric.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef Coord&lt; [**X**](structX.md), [**Y**](structY.md) &gt; | [**CoordArg**](#typedef-coordarg)  <br>_The type of a coordinate in the Cartesian coordinate system._  |
| typedef Coord&lt; Corner1Tag, Corner2Tag, Corner3Tag &gt; | [**CoordResult**](#typedef-coordresult)  <br>_The type of a coordinate in the barycentric coordinate system._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**CartesianToBarycentric**](#function-cartesiantobarycentric-13) (CartesianCoord const & corner1, CartesianCoord const & corner2, CartesianCoord const & corner3) <br>_Construct the operator which converts between the coordinate systems._  |
|   | [**CartesianToBarycentric**](#function-cartesiantobarycentric-23) ([**CartesianToBarycentric**](classCartesianToBarycentric.md) const & other) = default<br>_A copy operator for the mapping operator._  |
|   | [**CartesianToBarycentric**](#function-cartesiantobarycentric-33) ([**CartesianToBarycentric**](classCartesianToBarycentric.md) && x) = default<br>_A r-value copy operator for the mapping operator._  |
|  KOKKOS\_FUNCTION [**CoordResult**](classCartesianToBarycentric.md#typedef-coordresult) | [**operator()**](#function-operator) ([**CoordArg**](classCartesianToBarycentric.md#typedef-coordarg) const & pos) const<br>_The operator to get the equivalent barycentric coordinate of the Cartesian coordinate._  |
|  [**CartesianToBarycentric**](classCartesianToBarycentric.md) & | [**operator=**](#function-operator_1) ([**CartesianToBarycentric**](classCartesianToBarycentric.md) const & x) = default<br>_A copy operator for the mapping operator._  |
|  [**CartesianToBarycentric**](classCartesianToBarycentric.md) & | [**operator=**](#function-operator_2) ([**CartesianToBarycentric**](classCartesianToBarycentric.md) && x) = default<br>_A r-value copy operator for the mapping operator._  |
|   | [**~CartesianToBarycentric**](#function-cartesiantobarycentric) () = default<br>_The destructor of the mapping operator._  |




























## Detailed Description


Tags are used to identify the corners of the triangle. This ensures that there are different types for coordinate systems related to different triangles.




**Template parameters:**


* [**X**](structX.md) The tag of the x dimension of the Cartesian coordinates. 
* [**Y**](structY.md) The tag of the y dimension of the Cartesian coordinates. 
* `Corner1Tag` A tag identifying the first corner of the triangle. 
* `Corner2Tag` A tag identifying the second corner of the triangle. 
* `Corner3Tag` A tag identifying the third corner of the triangle. 




    
## Public Types Documentation




### typedef CoordArg 

_The type of a coordinate in the Cartesian coordinate system._ 
```C++
using CartesianToBarycentric< X, Y, Corner1Tag, Corner2Tag, Corner3Tag >::CoordArg =  Coord<X, Y>;
```




<hr>



### typedef CoordResult 

_The type of a coordinate in the barycentric coordinate system._ 
```C++
using CartesianToBarycentric< X, Y, Corner1Tag, Corner2Tag, Corner3Tag >::CoordResult =  Coord<Corner1Tag, Corner2Tag, Corner3Tag>;
```




<hr>
## Public Functions Documentation




### function CartesianToBarycentric [1/3]

_Construct the operator which converts between the coordinate systems._ 
```C++
inline CartesianToBarycentric::CartesianToBarycentric (
    CartesianCoord const & corner1,
    CartesianCoord const & corner2,
    CartesianCoord const & corner3
) 
```





**Parameters:**


* `corner1` The coordinates of the first corner of the triangle. 
* `corner2` The coordinates of the second corner of the triangle. 
* `corner3` The coordinates of the third corner of the triangle. 




        

<hr>



### function CartesianToBarycentric [2/3]

_A copy operator for the mapping operator._ 
```C++
CartesianToBarycentric::CartesianToBarycentric (
    CartesianToBarycentric const & other
) = default
```





**Parameters:**


* `other` The object to be copied. 




        

<hr>



### function CartesianToBarycentric [3/3]

_A r-value copy operator for the mapping operator._ 
```C++
CartesianToBarycentric::CartesianToBarycentric (
    CartesianToBarycentric && x
) = default
```





**Parameters:**


* `x` The object to be consumed. 




        

<hr>



### function operator() 

_The operator to get the equivalent barycentric coordinate of the Cartesian coordinate._ 
```C++
inline KOKKOS_FUNCTION CoordResult CartesianToBarycentric::operator() (
    CoordArg const & pos
) const
```





**Parameters:**


* `pos` The known Cartesian coordinate.



**Returns:**

The equivalent barycentric coordinate. 





        

<hr>



### function operator= 

_A copy operator for the mapping operator._ 
```C++
CartesianToBarycentric & CartesianToBarycentric::operator= (
    CartesianToBarycentric const & x
) = default
```





**Parameters:**


* `x` The object to be copied. 



**Returns:**

A reference to this class instance. 





        

<hr>



### function operator= 

_A r-value copy operator for the mapping operator._ 
```C++
CartesianToBarycentric & CartesianToBarycentric::operator= (
    CartesianToBarycentric && x
) = default
```





**Parameters:**


* `x` The object to be consumed. 



**Returns:**

A reference to this class instance. 





        

<hr>



### function ~CartesianToBarycentric 

_The destructor of the mapping operator._ 
```C++
CartesianToBarycentric::~CartesianToBarycentric () = default
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/mapping/cartesian_to_barycentric.hpp`

