

# Class CartesianToCylindrical

**template &lt;class [**X**](structX.md), class [**Y**](structY.md), class Z, class [**R**](structR.md), class Zeta&gt;**



[**ClassList**](annotated.md) **>** [**CartesianToCylindrical**](classCartesianToCylindrical.md)



_A class for describing the cylindrical 3D mapping._ [More...](#detailed-description)

* `#include <cartesian_to_cylindrical.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef Coord&lt; [**X**](structX.md), [**Y**](structY.md), Z &gt; | [**CoordArg**](#typedef-coordarg)  <br>_The type of the argument of the function described by this mapping._  |
| typedef Coord&lt; [**R**](structR.md), Z, Zeta &gt; | [**CoordResult**](#typedef-coordresult)  <br>_The type of the result of the function described by this mapping._  |
| typedef typename [**R::Dual**](structR.md#typedef-dual) | [**R\_cov**](#typedef-r_cov)  <br>_The covariant form of the radial cylindrical coordinate._  |
| typedef typename [**X::Dual**](structX.md#typedef-dual-14) | [**X\_cov**](#typedef-x_cov)  <br>_The covariant form of the first Cartesian coordinate._  |
| typedef typename [**Y::Dual**](structY.md#typedef-dual-13) | [**Y\_cov**](#typedef-y_cov)  <br>_The covariant form of the second Cartesian coordinate._  |
| typedef typename Z::Dual | [**Z\_cov**](#typedef-z_cov)  <br>_The covariant form of the third Cartesian coordinate._  |
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
|   | [**CartesianToCylindrical**](#function-cartesiantocylindrical-13) () = default<br> |
|  KOKKOS\_FUNCTION | [**CartesianToCylindrical**](#function-cartesiantocylindrical-23) ([**CartesianToCylindrical**](classCartesianToCylindrical.md) const & other) <br>_Instantiate a_ [_**CartesianToCylindrical**_](classCartesianToCylindrical.md) _from another_[_**CartesianToCylindrical**_](classCartesianToCylindrical.md) _(lvalue)._ |
|   | [**CartesianToCylindrical**](#function-cartesiantocylindrical-33) ([**CartesianToCylindrical**](classCartesianToCylindrical.md) && x) = default<br>_Instantiate a Curvilinear2DToCartesian from another temporary_ [_**CartesianToCylindrical**_](classCartesianToCylindrical.md) _(rvalue)._ |
|  [**CylindricalToCartesian**](classCylindricalToCartesian.md)&lt; [**R**](structR.md), Z, Zeta, [**X**](structX.md), [**Y**](structY.md) &gt; | [**get\_inverse\_mapping**](#function-get_inverse_mapping) () const<br>_Get the inverse mapping._  |
|  KOKKOS\_FUNCTION double | [**jacobian**](#function-jacobian) ([**CoordArg**](classCartesianToCylindrical.md#typedef-coordarg) const & coord) <br>_Compute the Jacobian, the determinant of the Jacobian matrix of the mapping._  |
|  KOKKOS\_FUNCTION double | [**jacobian\_component**](#function-jacobian_component) ([**CoordArg**](classCartesianToCylindrical.md#typedef-coordarg) const & coord) const<br>_Compute the (i,j) coefficient of the Jacobian matrix._  |
|  KOKKOS\_FUNCTION [**DTensor**](classTensor.md)&lt; VectorIndexSet&lt; [**R**](structR.md), Z, Zeta &gt;, VectorIndexSet&lt; [**X\_cov**](classCartesianToCylindrical.md#typedef-x_cov), [**Y\_cov**](classCartesianToCylindrical.md#typedef-y_cov), [**Z\_cov**](classCartesianToCylindrical.md#typedef-z_cov) &gt; &gt; | [**jacobian\_matrix**](#function-jacobian_matrix) ([**CoordArg**](classCartesianToCylindrical.md#typedef-coordarg) const & coord) const<br>_Compute full Jacobian matrix._  |
|  KOKKOS\_FUNCTION [**CoordResult**](classCartesianToCylindrical.md#typedef-coordresult) | [**operator()**](#function-operator) ([**CoordArg**](classCartesianToCylindrical.md#typedef-coordarg) const & coord) const<br>_Convert the coordinate (x,y,z) to the equivalent_  _coordinate._ |
|  [**CartesianToCylindrical**](classCartesianToCylindrical.md) & | [**operator=**](#function-operator_1) ([**CartesianToCylindrical**](classCartesianToCylindrical.md) const & x) = default<br>_Assign a_ [_**CartesianToCylindrical**_](classCartesianToCylindrical.md) _from another_[_**CartesianToCylindrical**_](classCartesianToCylindrical.md) _(lvalue)._ |
|  [**CartesianToCylindrical**](classCartesianToCylindrical.md) & | [**operator=**](#function-operator_2) ([**CartesianToCylindrical**](classCartesianToCylindrical.md) && x) = default<br>_Assign a_ [_**CartesianToCylindrical**_](classCartesianToCylindrical.md) _from another temporary_[_**CartesianToCylindrical**_](classCartesianToCylindrical.md) _(rvalue)._ |
|   | [**~CartesianToCylindrical**](#function-cartesiantocylindrical) () = default<br> |




























## Detailed Description


The mapping  is defined as follow :











It and its Jacobian matrix are invertible everywhere except for .


The Jacobian matrix coefficients are defined as follow





























and the matrix determinant: . 


    
## Public Types Documentation




### typedef CoordArg 

_The type of the argument of the function described by this mapping._ 
```C++
using CartesianToCylindrical< X, Y, Z, R, Zeta >::CoordArg =  Coord<X, Y, Z>;
```




<hr>



### typedef CoordResult 

_The type of the result of the function described by this mapping._ 
```C++
using CartesianToCylindrical< X, Y, Z, R, Zeta >::CoordResult =  Coord<R, Z, Zeta>;
```




<hr>



### typedef R\_cov 

_The covariant form of the radial cylindrical coordinate._ 
```C++
using CartesianToCylindrical< X, Y, Z, R, Zeta >::R_cov =  typename R::Dual;
```




<hr>



### typedef X\_cov 

_The covariant form of the first Cartesian coordinate._ 
```C++
using CartesianToCylindrical< X, Y, Z, R, Zeta >::X_cov =  typename X::Dual;
```




<hr>



### typedef Y\_cov 

_The covariant form of the second Cartesian coordinate._ 
```C++
using CartesianToCylindrical< X, Y, Z, R, Zeta >::Y_cov =  typename Y::Dual;
```




<hr>



### typedef Z\_cov 

_The covariant form of the third Cartesian coordinate._ 
```C++
using CartesianToCylindrical< X, Y, Z, R, Zeta >::Z_cov =  typename Z::Dual;
```




<hr>



### typedef Zeta\_cov 

_The covariant form of the angular cylindrical coordinate._ 
```C++
using CartesianToCylindrical< X, Y, Z, R, Zeta >::Zeta_cov =  typename Zeta::Dual;
```




<hr>



### typedef cartesian\_tag\_x 

_Indicate the first Cartesian coordinate._ 
```C++
using CartesianToCylindrical< X, Y, Z, R, Zeta >::cartesian_tag_x =  X;
```




<hr>



### typedef cartesian\_tag\_y 

_Indicate the second Cartesian coordinate._ 
```C++
using CartesianToCylindrical< X, Y, Z, R, Zeta >::cartesian_tag_y =  Y;
```




<hr>



### typedef cartesian\_tag\_z 

_Indicate the second Cartesian coordinate._ 
```C++
using CartesianToCylindrical< X, Y, Z, R, Zeta >::cartesian_tag_z =  Z;
```




<hr>



### typedef cylindrical\_tag\_R 

_Indicate the radial cylindrical coordinate._ 
```C++
using CartesianToCylindrical< X, Y, Z, R, Zeta >::cylindrical_tag_R =  R;
```




<hr>



### typedef cylindrical\_tag\_Z 

_Indicate the longitudinal cylindrical coordinate._ 
```C++
using CartesianToCylindrical< X, Y, Z, R, Zeta >::cylindrical_tag_Z =  Z;
```




<hr>



### typedef cylindrical\_tag\_Zeta 

_Indicate the angular cylindrical coordinate._ 
```C++
using CartesianToCylindrical< X, Y, Z, R, Zeta >::cylindrical_tag_Zeta =  Zeta;
```




<hr>
## Public Functions Documentation




### function CartesianToCylindrical [1/3]

```C++
CartesianToCylindrical::CartesianToCylindrical () = default
```




<hr>



### function CartesianToCylindrical [2/3]

_Instantiate a_ [_**CartesianToCylindrical**_](classCartesianToCylindrical.md) _from another_[_**CartesianToCylindrical**_](classCartesianToCylindrical.md) _(lvalue)._
```C++
inline KOKKOS_FUNCTION CartesianToCylindrical::CartesianToCylindrical (
    CartesianToCylindrical const & other
) 
```





**Parameters:**


* `other` [**CartesianToCylindrical**](classCartesianToCylindrical.md) mapping used to instantiate the new one. 




        

<hr>



### function CartesianToCylindrical [3/3]

_Instantiate a Curvilinear2DToCartesian from another temporary_ [_**CartesianToCylindrical**_](classCartesianToCylindrical.md) _(rvalue)._
```C++
CartesianToCylindrical::CartesianToCylindrical (
    CartesianToCylindrical && x
) = default
```





**Parameters:**


* `x` Curvilinear2DToCartesian mapping used to instantiate the new one. 




        

<hr>



### function get\_inverse\_mapping 

_Get the inverse mapping._ 
```C++
inline CylindricalToCartesian < R , Z, Zeta, X , Y > CartesianToCylindrical::get_inverse_mapping () const
```





**Returns:**

The inverse mapping. 





        

<hr>



### function jacobian 

_Compute the Jacobian, the determinant of the Jacobian matrix of the mapping._ 
```C++
inline KOKKOS_FUNCTION double CartesianToCylindrical::jacobian (
    CoordArg const & coord
) 
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
inline KOKKOS_FUNCTION double CartesianToCylindrical::jacobian_component (
    CoordArg const & coord
) const
```



Be careful because not all mappings are invertible, especially at the centre point.




**Parameters:**


* `coord` The coordinate where we evaluate the inverse Jacobian matrix.



**Returns:**

A double with the value of the (i,j) coefficient of the inverse Jacobian matrix. 





        

<hr>



### function jacobian\_matrix 

_Compute full Jacobian matrix._ 
```C++
inline KOKKOS_FUNCTION DTensor < VectorIndexSet< R , Z, Zeta >, VectorIndexSet< X_cov , Y_cov , Z_cov > > CartesianToCylindrical::jacobian_matrix (
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

_Convert the coordinate (x,y,z) to the equivalent_  _coordinate._
```C++
inline KOKKOS_FUNCTION CoordResult CartesianToCylindrical::operator() (
    CoordArg const & coord
) const
```





**Parameters:**


* `coord` The coordinate to be converted.



**Returns:**

The equivalent coordinate. 





        

<hr>



### function operator= 

_Assign a_ [_**CartesianToCylindrical**_](classCartesianToCylindrical.md) _from another_[_**CartesianToCylindrical**_](classCartesianToCylindrical.md) _(lvalue)._
```C++
CartesianToCylindrical & CartesianToCylindrical::operator= (
    CartesianToCylindrical const & x
) = default
```





**Parameters:**


* `x` [**CartesianToCylindrical**](classCartesianToCylindrical.md) mapping used to assign.



**Returns:**

The [**CartesianToCylindrical**](classCartesianToCylindrical.md) assigned. 





        

<hr>



### function operator= 

_Assign a_ [_**CartesianToCylindrical**_](classCartesianToCylindrical.md) _from another temporary_[_**CartesianToCylindrical**_](classCartesianToCylindrical.md) _(rvalue)._
```C++
CartesianToCylindrical & CartesianToCylindrical::operator= (
    CartesianToCylindrical && x
) = default
```





**Parameters:**


* `x` [**CartesianToCylindrical**](classCartesianToCylindrical.md) mapping used to assign.



**Returns:**

The [**CartesianToCylindrical**](classCartesianToCylindrical.md) assigned. 





        

<hr>



### function ~CartesianToCylindrical 

```C++
CartesianToCylindrical::~CartesianToCylindrical () = default
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/mapping/cartesian_to_cylindrical.hpp`

