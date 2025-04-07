

# Class BarycentricToCartesian

**template &lt;class Corner1Tag, class Corner2Tag, class Corner3Tag, class [**X**](structX.md), class [**Y**](structY.md)&gt;**



[**ClassList**](annotated.md) **>** [**BarycentricToCartesian**](classBarycentricToCartesian.md)



_A class to convert barycentric coordinates to Cartesian coordinates on a triangle._ [More...](#detailed-description)

* `#include <barycentric_to_cartesian.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef Coord&lt; Corner1Tag, Corner2Tag, Corner3Tag &gt; | [**CoordArg**](#typedef-coordarg)  <br>_The type of a coordinate in the barycentric coordinate system._  |
| typedef Coord&lt; [**X**](structX.md), [**Y**](structY.md) &gt; | [**CoordResult**](#typedef-coordresult)  <br>_The type of a coordinate in the Cartesian coordinate system._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**BarycentricToCartesian**](#function-barycentrictocartesian-13) (CartesianCoord const & corner1, CartesianCoord const & corner2, CartesianCoord const & corner3) <br>_Construct the operator which converts between the coordinate systems._  |
|   | [**BarycentricToCartesian**](#function-barycentrictocartesian-23) ([**BarycentricToCartesian**](classBarycentricToCartesian.md) const & other) = default<br>_A copy operator for the mapping operator._  |
|   | [**BarycentricToCartesian**](#function-barycentrictocartesian-33) ([**BarycentricToCartesian**](classBarycentricToCartesian.md) && x) = default<br>_A r-value copy operator for the mapping operator._  |
|  [**CartesianToBarycentric**](classCartesianToBarycentric.md)&lt; [**X**](structX.md), [**Y**](structY.md), Corner1Tag, Corner2Tag, Corner3Tag &gt; | [**get\_inverse\_mapping**](#function-get_inverse_mapping) () const<br>_Get the inverse mapping._  |
|  KOKKOS\_FUNCTION [**CoordResult**](classBarycentricToCartesian.md#typedef-coordresult) | [**operator()**](#function-operator) ([**CoordArg**](classBarycentricToCartesian.md#typedef-coordarg) const & pos) const<br>_The operator to get the equivalent Cartesian coordinate of the barycentric coordinate._  |
|  [**BarycentricToCartesian**](classBarycentricToCartesian.md) & | [**operator=**](#function-operator_1) ([**BarycentricToCartesian**](classBarycentricToCartesian.md) const & x) = default<br>_A copy operator for the mapping operator._  |
|  [**BarycentricToCartesian**](classBarycentricToCartesian.md) & | [**operator=**](#function-operator_2) ([**BarycentricToCartesian**](classBarycentricToCartesian.md) && x) = default<br>_A r-value copy operator for the mapping operator._  |
|   | [**~BarycentricToCartesian**](#function-barycentrictocartesian) () = default<br>_The destructor of the mapping operator._  |




























## Detailed Description


Tags are used to identify the corners of the triangle. This ensures that there are different types for coordinate systems related to different triangles.




**Template parameters:**


* `Corner1Tag` A tag identifying the first corner of the triangle. 
* `Corner2Tag` A tag identifying the second corner of the triangle. 
* `Corner3Tag` A tag identifying the third corner of the triangle. 
* [**X**](structX.md) The tag of the x dimension of the Cartesian coordinates. 
* [**Y**](structY.md) The tag of the y dimension of the Cartesian coordinates. 




    
## Public Types Documentation




### typedef CoordArg 

_The type of a coordinate in the barycentric coordinate system._ 
```C++
using BarycentricToCartesian< Corner1Tag, Corner2Tag, Corner3Tag, X, Y >::CoordArg =  Coord<Corner1Tag, Corner2Tag, Corner3Tag>;
```




<hr>



### typedef CoordResult 

_The type of a coordinate in the Cartesian coordinate system._ 
```C++
using BarycentricToCartesian< Corner1Tag, Corner2Tag, Corner3Tag, X, Y >::CoordResult =  Coord<X, Y>;
```




<hr>
## Public Functions Documentation




### function BarycentricToCartesian [1/3]

_Construct the operator which converts between the coordinate systems._ 
```C++
inline BarycentricToCartesian::BarycentricToCartesian (
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



### function BarycentricToCartesian [2/3]

_A copy operator for the mapping operator._ 
```C++
BarycentricToCartesian::BarycentricToCartesian (
    BarycentricToCartesian const & other
) = default
```





**Parameters:**


* `other` The object to be copied. 




        

<hr>



### function BarycentricToCartesian [3/3]

_A r-value copy operator for the mapping operator._ 
```C++
BarycentricToCartesian::BarycentricToCartesian (
    BarycentricToCartesian && x
) = default
```





**Parameters:**


* `x` The object to be consumed. 




        

<hr>



### function get\_inverse\_mapping 

_Get the inverse mapping._ 
```C++
inline CartesianToBarycentric < X , Y , Corner1Tag, Corner2Tag, Corner3Tag > BarycentricToCartesian::get_inverse_mapping () const
```





**Returns:**

The inverse mapping. 





        

<hr>



### function operator() 

_The operator to get the equivalent Cartesian coordinate of the barycentric coordinate._ 
```C++
inline KOKKOS_FUNCTION CoordResult BarycentricToCartesian::operator() (
    CoordArg const & pos
) const
```





**Parameters:**


* `pos` The known barycentric coordinate.



**Returns:**

The equivalent Cartesian coordinate. 





        

<hr>



### function operator= 

_A copy operator for the mapping operator._ 
```C++
BarycentricToCartesian & BarycentricToCartesian::operator= (
    BarycentricToCartesian const & x
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
BarycentricToCartesian & BarycentricToCartesian::operator= (
    BarycentricToCartesian && x
) = default
```





**Parameters:**


* `x` The object to be consumed. 



**Returns:**

A reference to this class instance. 





        

<hr>



### function ~BarycentricToCartesian 

_The destructor of the mapping operator._ 
```C++
BarycentricToCartesian::~BarycentricToCartesian () = default
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/mapping/barycentric_to_cartesian.hpp`

