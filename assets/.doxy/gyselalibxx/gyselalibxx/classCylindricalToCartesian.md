

# Class CylindricalToCartesian

**template &lt;class [**R**](structR.md), class Z, class Zeta, class [**X**](structX.md), class [**Y**](structY.md)&gt;**



[**ClassList**](annotated.md) **>** [**CylindricalToCartesian**](classCylindricalToCartesian.md)



_A class for describing the cylindrical 3D mapping._ [More...](#detailed-description)

* `#include <cylindrical_to_cartesian.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef Coord&lt; [**R**](structR.md), Z, Zeta &gt; | [**CoordArg**](#typedef-coordarg)  <br>_The type of the argument of the function described by this mapping._  |
| typedef Coord&lt; [**X**](structX.md), [**Y**](structY.md), Z &gt; | [**CoordResult**](#typedef-coordresult)  <br>_The type of the result of the function described by this mapping._  |
| typedef typename [**R::Dual**](structR.md#typedef-dual) | [**R\_cov**](#typedef-r_cov)  <br>_The covariant form of the radial cylindrical coordinate._  |
| typedef typename [**X::Dual**](structX.md#typedef-dual-14) | [**X\_cov**](#typedef-x_cov)  <br>_The covariant form of the first Cartesian coordinate._  |
| typedef typename [**Y::Dual**](structY.md#typedef-dual-13) | [**Y\_cov**](#typedef-y_cov)  <br>_The covariant form of the second Cartesian coordinate._  |
| typedef typename Z::Dual | [**Z\_cov**](#typedef-z_cov)  <br>_The covariant form of the third Cartesian coordinate and the longitudinal cylindrical coordinate._  |
| typedef typename Zeta::Dual | [**Zeta\_cov**](#typedef-zeta_cov)  <br>_The covariant form of the angular cylindrical coordinate._  |
| typedef [**X**](structX.md) | [**cartesian\_tag\_x**](#typedef-cartesian_tag_x)  <br>_Indicate the first Cartesian coordinate._  |
| typedef [**Y**](structY.md) | [**cartesian\_tag\_y**](#typedef-cartesian_tag_y)  <br>_Indicate the second Cartesian coordinate._  |
| typedef Z | [**cartesian\_tag\_z**](#typedef-cartesian_tag_z)  <br>_Indicate the second Cartesian coordinate._  |
| typedef [**R**](structR.md) | [**cylindrical\_tag\_R**](#typedef-cylindrical_tag_r)  <br>_Indicate the radial cylindrical coordinate._  |
| typedef Z | [**cylindrical\_tag\_Z**](#typedef-cylindrical_tag_z)  <br>_Indicate the longitudinal cylindrical coordinate._  |
| typedef Zeta | [**cylindrical\_tag\_Zeta**](#typedef-cylindrical_tag_zeta)  <br>_Indicate the angular cylindrical coordinate._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**CylindricalToCartesian**](#function-cylindricaltocartesian-13) () = default<br> |
|  KOKKOS\_FUNCTION | [**CylindricalToCartesian**](#function-cylindricaltocartesian-23) ([**CylindricalToCartesian**](classCylindricalToCartesian.md) const & other) <br>_Instantiate a_ [_**CylindricalToCartesian**_](classCylindricalToCartesian.md) _from another_[_**CylindricalToCartesian**_](classCylindricalToCartesian.md) _(lvalue)._ |
|   | [**CylindricalToCartesian**](#function-cylindricaltocartesian-33) ([**CylindricalToCartesian**](classCylindricalToCartesian.md) && x) = default<br>_Instantiate a_ [_**CylindricalToCartesian**_](classCylindricalToCartesian.md) _from another temporary_[_**CylindricalToCartesian**_](classCylindricalToCartesian.md) _(rvalue)._ |
|  [**CartesianToCylindrical**](classCartesianToCylindrical.md)&lt; [**X**](structX.md), [**Y**](structY.md), Z, [**R**](structR.md), Zeta &gt; | [**get\_inverse\_mapping**](#function-get_inverse_mapping) () const<br>_Get the inverse mapping._  |
|  KOKKOS\_FUNCTION double | [**inv\_jacobian\_component**](#function-inv_jacobian_component) ([**CoordArg**](classCylindricalToCartesian.md#typedef-coordarg) const & coord) const<br>_Compute the (i,j) coefficient of the inverse Jacobian matrix._  |
|  KOKKOS\_FUNCTION [**DTensor**](classTensor.md)&lt; VectorIndexSet&lt; [**R**](structR.md), Z, Zeta &gt;, VectorIndexSet&lt; [**X\_cov**](classCylindricalToCartesian.md#typedef-x_cov), [**Y\_cov**](classCylindricalToCartesian.md#typedef-y_cov), [**Z\_cov**](classCylindricalToCartesian.md#typedef-z_cov) &gt; &gt; | [**inv\_jacobian\_matrix**](#function-inv_jacobian_matrix) ([**CoordArg**](classCylindricalToCartesian.md#typedef-coordarg) const & coord) const<br>_Compute full inverse Jacobian matrix._  |
|  KOKKOS\_FUNCTION double | [**jacobian**](#function-jacobian) ([**CoordArg**](classCylindricalToCartesian.md#typedef-coordarg) const & coord) const<br>_Compute the Jacobian, the determinant of the Jacobian matrix of the mapping._  |
|  KOKKOS\_FUNCTION double | [**jacobian\_component**](#function-jacobian_component) ([**CoordArg**](classCylindricalToCartesian.md#typedef-coordarg) const & coord) const<br>_Compute the (i,j) coefficient of the Jacobian matrix._  |
|  KOKKOS\_FUNCTION [**DTensor**](classTensor.md)&lt; VectorIndexSet&lt; [**X**](structX.md), [**Y**](structY.md), Z &gt;, VectorIndexSet&lt; [**R\_cov**](classCylindricalToCartesian.md#typedef-r_cov), [**Z\_cov**](classCylindricalToCartesian.md#typedef-z_cov), [**Zeta\_cov**](classCylindricalToCartesian.md#typedef-zeta_cov) &gt; &gt; | [**jacobian\_matrix**](#function-jacobian_matrix) ([**CoordArg**](classCylindricalToCartesian.md#typedef-coordarg) const & coord) const<br>_Compute full Jacobian matrix._  |
|  KOKKOS\_FUNCTION [**CoordResult**](classCylindricalToCartesian.md#typedef-coordresult) | [**operator()**](#function-operator) ([**CoordArg**](classCylindricalToCartesian.md#typedef-coordarg) const & coord) const<br>_Convert the_  _coordinate to the equivalent (x,y) coordinate._ |
|  [**CylindricalToCartesian**](classCylindricalToCartesian.md) & | [**operator=**](#function-operator_1) ([**CylindricalToCartesian**](classCylindricalToCartesian.md) const & x) = default<br>_Assign a_ [_**CylindricalToCartesian**_](classCylindricalToCartesian.md) _from another_[_**CylindricalToCartesian**_](classCylindricalToCartesian.md) _(lvalue)._ |
|  [**CylindricalToCartesian**](classCylindricalToCartesian.md) & | [**operator=**](#function-operator_2) ([**CylindricalToCartesian**](classCylindricalToCartesian.md) && x) = default<br>_Assign a_ [_**CylindricalToCartesian**_](classCylindricalToCartesian.md) _from another temporary_[_**CylindricalToCartesian**_](classCylindricalToCartesian.md) _(rvalue)._ |
|   | [**~CylindricalToCartesian**](#function-cylindricaltocartesian) () = default<br> |




























## Detailed Description


The mapping  is defined as follow :











It and its Jacobian matrix are invertible everywhere except for .


The Jacobian matrix coefficients are defined as follow





























and the matrix determinant: . 


    
## Public Types Documentation




### typedef CoordArg 

_The type of the argument of the function described by this mapping._ 
```C++
using CylindricalToCartesian< R, Z, Zeta, X, Y >::CoordArg =  Coord<R, Z, Zeta>;
```




<hr>



### typedef CoordResult 

_The type of the result of the function described by this mapping._ 
```C++
using CylindricalToCartesian< R, Z, Zeta, X, Y >::CoordResult =  Coord<X, Y, Z>;
```




<hr>



### typedef R\_cov 

_The covariant form of the radial cylindrical coordinate._ 
```C++
using CylindricalToCartesian< R, Z, Zeta, X, Y >::R_cov =  typename R::Dual;
```




<hr>



### typedef X\_cov 

_The covariant form of the first Cartesian coordinate._ 
```C++
using CylindricalToCartesian< R, Z, Zeta, X, Y >::X_cov =  typename X::Dual;
```




<hr>



### typedef Y\_cov 

_The covariant form of the second Cartesian coordinate._ 
```C++
using CylindricalToCartesian< R, Z, Zeta, X, Y >::Y_cov =  typename Y::Dual;
```




<hr>



### typedef Z\_cov 

_The covariant form of the third Cartesian coordinate and the longitudinal cylindrical coordinate._ 
```C++
using CylindricalToCartesian< R, Z, Zeta, X, Y >::Z_cov =  typename Z::Dual;
```




<hr>



### typedef Zeta\_cov 

_The covariant form of the angular cylindrical coordinate._ 
```C++
using CylindricalToCartesian< R, Z, Zeta, X, Y >::Zeta_cov =  typename Zeta::Dual;
```




<hr>



### typedef cartesian\_tag\_x 

_Indicate the first Cartesian coordinate._ 
```C++
using CylindricalToCartesian< R, Z, Zeta, X, Y >::cartesian_tag_x =  X;
```




<hr>



### typedef cartesian\_tag\_y 

_Indicate the second Cartesian coordinate._ 
```C++
using CylindricalToCartesian< R, Z, Zeta, X, Y >::cartesian_tag_y =  Y;
```




<hr>



### typedef cartesian\_tag\_z 

_Indicate the second Cartesian coordinate._ 
```C++
using CylindricalToCartesian< R, Z, Zeta, X, Y >::cartesian_tag_z =  Z;
```




<hr>



### typedef cylindrical\_tag\_R 

_Indicate the radial cylindrical coordinate._ 
```C++
using CylindricalToCartesian< R, Z, Zeta, X, Y >::cylindrical_tag_R =  R;
```




<hr>



### typedef cylindrical\_tag\_Z 

_Indicate the longitudinal cylindrical coordinate._ 
```C++
using CylindricalToCartesian< R, Z, Zeta, X, Y >::cylindrical_tag_Z =  Z;
```




<hr>



### typedef cylindrical\_tag\_Zeta 

_Indicate the angular cylindrical coordinate._ 
```C++
using CylindricalToCartesian< R, Z, Zeta, X, Y >::cylindrical_tag_Zeta =  Zeta;
```




<hr>
## Public Functions Documentation




### function CylindricalToCartesian [1/3]

```C++
CylindricalToCartesian::CylindricalToCartesian () = default
```




<hr>



### function CylindricalToCartesian [2/3]

_Instantiate a_ [_**CylindricalToCartesian**_](classCylindricalToCartesian.md) _from another_[_**CylindricalToCartesian**_](classCylindricalToCartesian.md) _(lvalue)._
```C++
inline KOKKOS_FUNCTION CylindricalToCartesian::CylindricalToCartesian (
    CylindricalToCartesian const & other
) 
```





**Parameters:**


* `other` [**CylindricalToCartesian**](classCylindricalToCartesian.md) mapping used to instantiate the new one. 




        

<hr>



### function CylindricalToCartesian [3/3]

_Instantiate a_ [_**CylindricalToCartesian**_](classCylindricalToCartesian.md) _from another temporary_[_**CylindricalToCartesian**_](classCylindricalToCartesian.md) _(rvalue)._
```C++
CylindricalToCartesian::CylindricalToCartesian (
    CylindricalToCartesian && x
) = default
```





**Parameters:**


* `x` [**CylindricalToCartesian**](classCylindricalToCartesian.md) mapping used to instantiate the new one. 




        

<hr>



### function get\_inverse\_mapping 

_Get the inverse mapping._ 
```C++
inline CartesianToCylindrical < X , Y , Z, R , Zeta > CylindricalToCartesian::get_inverse_mapping () const
```





**Returns:**

The inverse mapping. 





        

<hr>



### function inv\_jacobian\_component 

_Compute the (i,j) coefficient of the inverse Jacobian matrix._ 
```C++
template<class IndexTag1, class IndexTag2>
inline KOKKOS_FUNCTION double CylindricalToCartesian::inv_jacobian_component (
    CoordArg const & coord
) const
```



Be careful because not all mappings are invertible, especially at the centre point.




**Parameters:**


* `coord` The coordinate where we evaluate the inverse Jacobian matrix.



**Returns:**

A double with the value of the (i,j) coefficient of the inverse Jacobian matrix. 





        

<hr>



### function inv\_jacobian\_matrix 

_Compute full inverse Jacobian matrix._ 
```C++
inline KOKKOS_FUNCTION DTensor < VectorIndexSet< R , Z, Zeta >, VectorIndexSet< X_cov , Y_cov , Z_cov > > CylindricalToCartesian::inv_jacobian_matrix (
    CoordArg const & coord
) const
```



For some computations, we need the complete inverse Jacobian matrix or just the coefficients. The coefficients can be given independently with the function inv\_jacobian\_component.




**Parameters:**


* `coord` The coordinate where we evaluate the Jacobian matrix. 



**Returns:**

The inverse Jacobian matrix. 





        

<hr>



### function jacobian 

_Compute the Jacobian, the determinant of the Jacobian matrix of the mapping._ 
```C++
inline KOKKOS_FUNCTION double CylindricalToCartesian::jacobian (
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
inline KOKKOS_FUNCTION double CylindricalToCartesian::jacobian_component (
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
inline KOKKOS_FUNCTION DTensor < VectorIndexSet< X , Y , Z >, VectorIndexSet< R_cov , Z_cov , Zeta_cov > > CylindricalToCartesian::jacobian_matrix (
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

_Convert the_  _coordinate to the equivalent (x,y) coordinate._
```C++
inline KOKKOS_FUNCTION CoordResult CylindricalToCartesian::operator() (
    CoordArg const & coord
) const
```





**Parameters:**


* `coord` The coordinate to be converted.



**Returns:**

The equivalent coordinate. 





        

<hr>



### function operator= 

_Assign a_ [_**CylindricalToCartesian**_](classCylindricalToCartesian.md) _from another_[_**CylindricalToCartesian**_](classCylindricalToCartesian.md) _(lvalue)._
```C++
CylindricalToCartesian & CylindricalToCartesian::operator= (
    CylindricalToCartesian const & x
) = default
```





**Parameters:**


* `x` [**CylindricalToCartesian**](classCylindricalToCartesian.md) mapping used to assign.



**Returns:**

The [**CylindricalToCartesian**](classCylindricalToCartesian.md) assigned. 





        

<hr>



### function operator= 

_Assign a_ [_**CylindricalToCartesian**_](classCylindricalToCartesian.md) _from another temporary_[_**CylindricalToCartesian**_](classCylindricalToCartesian.md) _(rvalue)._
```C++
CylindricalToCartesian & CylindricalToCartesian::operator= (
    CylindricalToCartesian && x
) = default
```





**Parameters:**


* `x` [**CylindricalToCartesian**](classCylindricalToCartesian.md) mapping used to assign.



**Returns:**

The [**CylindricalToCartesian**](classCylindricalToCartesian.md) assigned. 





        

<hr>



### function ~CylindricalToCartesian 

```C++
CylindricalToCartesian::~CylindricalToCartesian () = default
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/mapping/cylindrical_to_cartesian.hpp`

