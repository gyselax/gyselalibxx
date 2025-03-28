

# Class CartesianToCircular

**template &lt;class [**X**](structX.md), class [**Y**](structY.md), class [**R**](structR.md), class [**Theta**](structTheta.md)&gt;**



[**ClassList**](annotated.md) **>** [**CartesianToCircular**](classCartesianToCircular.md)



_A class for describing the circular 2D mapping._ [More...](#detailed-description)

* `#include <cartesian_to_circular.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef Coord&lt; [**X**](structX.md), [**Y**](structY.md) &gt; | [**CoordArg**](#typedef-coordarg)  <br>_The type of the argument of the function described by this mapping._  |
| typedef Coord&lt; [**R**](structR.md), [**Theta**](structTheta.md) &gt; | [**CoordResult**](#typedef-coordresult)  <br>_The type of the result of the function described by this mapping._  |
| typedef typename [**R::Dual**](structR.md#typedef-dual) | [**R\_cov**](#typedef-r_cov)  <br>_The covariant form of the first logical coordinate._  |
| typedef typename [**Theta::Dual**](structTheta.md#typedef-dual) | [**Theta\_cov**](#typedef-theta_cov)  <br>_The covariant form of the second logical coordinate._  |
| typedef typename [**X::Dual**](structX.md#typedef-dual-14) | [**X\_cov**](#typedef-x_cov)  <br>_The covariant form of the first physical coordinate._  |
| typedef typename [**Y::Dual**](structY.md#typedef-dual-13) | [**Y\_cov**](#typedef-y_cov)  <br>_The covariant form of the second physical coordinate._  |
| typedef [**X**](structX.md) | [**cartesian\_tag\_x**](#typedef-cartesian_tag_x)  <br>_Indicate the first physical coordinate._  |
| typedef [**Y**](structY.md) | [**cartesian\_tag\_y**](#typedef-cartesian_tag_y)  <br>_Indicate the second physical coordinate._  |
| typedef [**R**](structR.md) | [**curvilinear\_tag\_r**](#typedef-curvilinear_tag_r)  <br>_Indicate the first logical coordinate._  |
| typedef [**Theta**](structTheta.md) | [**curvilinear\_tag\_theta**](#typedef-curvilinear_tag_theta)  <br>_Indicate the second logical coordinate._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**CartesianToCircular**](#function-cartesiantocircular-13) () = default<br> |
|  KOKKOS\_FUNCTION | [**CartesianToCircular**](#function-cartesiantocircular-23) ([**CartesianToCircular**](classCartesianToCircular.md) const & other) <br>_Instantiate a_ [_**CartesianToCircular**_](classCartesianToCircular.md) _from another_[_**CartesianToCircular**_](classCartesianToCircular.md) _(lvalue)._ |
|   | [**CartesianToCircular**](#function-cartesiantocircular-33) ([**CartesianToCircular**](classCartesianToCircular.md) && x) = default<br>_Instantiate a Curvilinear2DToCartesian from another temporary_ [_**CartesianToCircular**_](classCartesianToCircular.md) _(rvalue)._ |
|  [**CircularToCartesian**](classCircularToCartesian.md)&lt; [**R**](structR.md), [**Theta**](structTheta.md), [**X**](structX.md), [**Y**](structY.md) &gt; | [**get\_inverse\_mapping**](#function-get_inverse_mapping) () const<br>_Get the inverse mapping._  |
|  KOKKOS\_FUNCTION double | [**jacobian**](#function-jacobian) (Coord&lt; [**X**](structX.md), [**Y**](structY.md) &gt; const & coord) <br>_Compute the Jacobian, the determinant of the Jacobian matrix of the mapping._  |
|  KOKKOS\_INLINE\_FUNCTION double | [**jacobian\_component**](#function-jacobian_component) (Coord&lt; [**X**](structX.md), [**Y**](structY.md) &gt; const & coord) const<br>_Compute the (i,j) coefficient of the Jacobian matrix._  |
|  KOKKOS\_FUNCTION [**DTensor**](classTensor.md)&lt; VectorIndexSet&lt; [**X**](structX.md), [**Y**](structY.md) &gt;, VectorIndexSet&lt; [**R\_cov**](classCartesianToCircular.md#typedef-r_cov), [**Theta\_cov**](classCartesianToCircular.md#typedef-theta_cov) &gt; &gt; | [**jacobian\_matrix**](#function-jacobian_matrix) (Coord&lt; [**X**](structX.md), [**Y**](structY.md) &gt; const & coord) const<br>_Compute full Jacobian matrix._  |
|  KOKKOS\_FUNCTION Coord&lt; [**R**](structR.md), [**Theta**](structTheta.md) &gt; | [**operator()**](#function-operator) (Coord&lt; [**X**](structX.md), [**Y**](structY.md) &gt; const & coord) const<br>_Convert the coordinate (x,y) to the equivalent_  _coordinate._ |
|  [**CartesianToCircular**](classCartesianToCircular.md) & | [**operator=**](#function-operator_1) ([**CartesianToCircular**](classCartesianToCircular.md) const & x) = default<br>_Assign a_ [_**CartesianToCircular**_](classCartesianToCircular.md) _from another_[_**CartesianToCircular**_](classCartesianToCircular.md) _(lvalue)._ |
|  [**CartesianToCircular**](classCartesianToCircular.md) & | [**operator=**](#function-operator_2) ([**CartesianToCircular**](classCartesianToCircular.md) && x) = default<br>_Assign a_ [_**CartesianToCircular**_](classCartesianToCircular.md) _from another temporary_[_**CartesianToCircular**_](classCartesianToCircular.md) _(rvalue)._ |
|   | [**~CartesianToCircular**](#function-cartesiantocircular) () = default<br> |




























## Detailed Description


The mapping  is defined as follow :








It and its Jacobian matrix are invertible everywhere except for .


The Jacobian matrix coefficients are defined as follow














and the matrix determinant: . 


    
## Public Types Documentation




### typedef CoordArg 

_The type of the argument of the function described by this mapping._ 
```C++
using CartesianToCircular< X, Y, R, Theta >::CoordArg =  Coord<X, Y>;
```




<hr>



### typedef CoordResult 

_The type of the result of the function described by this mapping._ 
```C++
using CartesianToCircular< X, Y, R, Theta >::CoordResult =  Coord<R, Theta>;
```




<hr>



### typedef R\_cov 

_The covariant form of the first logical coordinate._ 
```C++
using CartesianToCircular< X, Y, R, Theta >::R_cov =  typename R::Dual;
```




<hr>



### typedef Theta\_cov 

_The covariant form of the second logical coordinate._ 
```C++
using CartesianToCircular< X, Y, R, Theta >::Theta_cov =  typename Theta::Dual;
```




<hr>



### typedef X\_cov 

_The covariant form of the first physical coordinate._ 
```C++
using CartesianToCircular< X, Y, R, Theta >::X_cov =  typename X::Dual;
```




<hr>



### typedef Y\_cov 

_The covariant form of the second physical coordinate._ 
```C++
using CartesianToCircular< X, Y, R, Theta >::Y_cov =  typename Y::Dual;
```




<hr>



### typedef cartesian\_tag\_x 

_Indicate the first physical coordinate._ 
```C++
using CartesianToCircular< X, Y, R, Theta >::cartesian_tag_x =  X;
```




<hr>



### typedef cartesian\_tag\_y 

_Indicate the second physical coordinate._ 
```C++
using CartesianToCircular< X, Y, R, Theta >::cartesian_tag_y =  Y;
```




<hr>



### typedef curvilinear\_tag\_r 

_Indicate the first logical coordinate._ 
```C++
using CartesianToCircular< X, Y, R, Theta >::curvilinear_tag_r =  R;
```




<hr>



### typedef curvilinear\_tag\_theta 

_Indicate the second logical coordinate._ 
```C++
using CartesianToCircular< X, Y, R, Theta >::curvilinear_tag_theta =  Theta;
```




<hr>
## Public Functions Documentation




### function CartesianToCircular [1/3]

```C++
CartesianToCircular::CartesianToCircular () = default
```




<hr>



### function CartesianToCircular [2/3]

_Instantiate a_ [_**CartesianToCircular**_](classCartesianToCircular.md) _from another_[_**CartesianToCircular**_](classCartesianToCircular.md) _(lvalue)._
```C++
inline KOKKOS_FUNCTION CartesianToCircular::CartesianToCircular (
    CartesianToCircular const & other
) 
```





**Parameters:**


* `other` [**CartesianToCircular**](classCartesianToCircular.md) mapping used to instantiate the new one. 




        

<hr>



### function CartesianToCircular [3/3]

_Instantiate a Curvilinear2DToCartesian from another temporary_ [_**CartesianToCircular**_](classCartesianToCircular.md) _(rvalue)._
```C++
CartesianToCircular::CartesianToCircular (
    CartesianToCircular && x
) = default
```





**Parameters:**


* `x` Curvilinear2DToCartesian mapping used to instantiate the new one. 




        

<hr>



### function get\_inverse\_mapping 

_Get the inverse mapping._ 
```C++
inline CircularToCartesian < R , Theta , X , Y > CartesianToCircular::get_inverse_mapping () const
```





**Returns:**

The inverse mapping. 





        

<hr>



### function jacobian 

_Compute the Jacobian, the determinant of the Jacobian matrix of the mapping._ 
```C++
inline KOKKOS_FUNCTION double CartesianToCircular::jacobian (
    Coord< X , Y > const & coord
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
inline KOKKOS_INLINE_FUNCTION double CartesianToCircular::jacobian_component (
    Coord< X , Y > const & coord
) const
```





**Parameters:**


* `coord` The coordinate where we evaluate the Jacobian matrix.



**Returns:**

A double with the value of the (1,1) coefficient of the Jacobian matrix. 





        

<hr>



### function jacobian\_matrix 

_Compute full Jacobian matrix._ 
```C++
inline KOKKOS_FUNCTION DTensor < VectorIndexSet< X , Y >, VectorIndexSet< R_cov , Theta_cov > > CartesianToCircular::jacobian_matrix (
    Coord< X , Y > const & coord
) const
```



For some computations, we need the complete Jacobian matrix or just the coefficients. The coefficients can be given independently with the function jacobian\_component.




**Parameters:**


* `coord` The coordinate where we evaluate the Jacobian matrix. 



**Returns:**

The Jacobian matrix. 





        

<hr>



### function operator() 

_Convert the coordinate (x,y) to the equivalent_  _coordinate._
```C++
inline KOKKOS_FUNCTION Coord< R , Theta > CartesianToCircular::operator() (
    Coord< X , Y > const & coord
) const
```





**Parameters:**


* `coord` The coordinate to be converted.



**Returns:**

The equivalent coordinate. 





        

<hr>



### function operator= 

_Assign a_ [_**CartesianToCircular**_](classCartesianToCircular.md) _from another_[_**CartesianToCircular**_](classCartesianToCircular.md) _(lvalue)._
```C++
CartesianToCircular & CartesianToCircular::operator= (
    CartesianToCircular const & x
) = default
```





**Parameters:**


* `x` [**CartesianToCircular**](classCartesianToCircular.md) mapping used to assign.



**Returns:**

The [**CartesianToCircular**](classCartesianToCircular.md) assigned. 





        

<hr>



### function operator= 

_Assign a_ [_**CartesianToCircular**_](classCartesianToCircular.md) _from another temporary_[_**CartesianToCircular**_](classCartesianToCircular.md) _(rvalue)._
```C++
CartesianToCircular & CartesianToCircular::operator= (
    CartesianToCircular && x
) = default
```





**Parameters:**


* `x` [**CartesianToCircular**](classCartesianToCircular.md) mapping used to assign.



**Returns:**

The [**CartesianToCircular**](classCartesianToCircular.md) assigned. 





        

<hr>



### function ~CartesianToCircular 

```C++
CartesianToCircular::~CartesianToCircular () = default
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/mapping/cartesian_to_circular.hpp`

