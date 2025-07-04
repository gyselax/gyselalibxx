

# Class CircularToCartesian

**template &lt;class [**R**](structR.md), class [**Theta**](structTheta.md), class [**X**](structX.md), class [**Y**](structY.md)&gt;**



[**ClassList**](annotated.md) **>** [**CircularToCartesian**](classCircularToCartesian.md)



_A class for describing the circular 2D mapping._ [More...](#detailed-description)

* `#include <circular_to_cartesian.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef Coord&lt; [**R**](structR.md), [**Theta**](structTheta.md) &gt; | [**CoordArg**](#typedef-coordarg)  <br>_The type of the argument of the function described by this mapping._  |
| typedef [**CoordArg**](classCircularToCartesian.md#typedef-coordarg) | [**CoordJacobian**](#typedef-coordjacobian)  <br>_The type of the coordinate that can be used to evaluate the Jacobian of this mapping._  |
| typedef Coord&lt; [**X**](structX.md), [**Y**](structY.md) &gt; | [**CoordResult**](#typedef-coordresult)  <br>_The type of the result of the function described by this mapping._  |
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
|  KOKKOS\_FUNCTION | [**CircularToCartesian**](#function-circulartocartesian-13) (Coord&lt; [**X**](structX.md), [**Y**](structY.md) &gt; o\_point=Coord&lt; [**X**](structX.md), [**Y**](structY.md) &gt;(0.0, 0.0)) <br>_Instantiate a_ [_**CircularToCartesian**_](classCircularToCartesian.md) _from parameters._ |
|  KOKKOS\_DEFAULTED\_FUNCTION | [**CircularToCartesian**](#function-circulartocartesian-23) ([**CircularToCartesian**](classCircularToCartesian.md) const & other) = default<br>_Instantiate a_ [_**CircularToCartesian**_](classCircularToCartesian.md) _from another_[_**CircularToCartesian**_](classCircularToCartesian.md) _(lvalue)._ |
|   | [**CircularToCartesian**](#function-circulartocartesian-33) ([**CircularToCartesian**](classCircularToCartesian.md) && x) = default<br>_Instantiate a_ [_**CircularToCartesian**_](classCircularToCartesian.md) _from another temporary_[_**CircularToCartesian**_](classCircularToCartesian.md) _(rvalue)._ |
|  KOKKOS\_INLINE\_FUNCTION [**CartesianToCircular**](classCartesianToCircular.md)&lt; [**X**](structX.md), [**Y**](structY.md), [**R**](structR.md), [**Theta**](structTheta.md) &gt; | [**get\_inverse\_mapping**](#function-get_inverse_mapping) () const<br>_Get the inverse mapping._  |
|  KOKKOS\_FUNCTION double | [**inv\_jacobian\_component**](#function-inv_jacobian_component) (Coord&lt; [**R**](structR.md), [**Theta**](structTheta.md) &gt; const & coord) const<br>_Compute the (i,j) coefficient of the inverse Jacobian matrix._  |
|  KOKKOS\_FUNCTION [**DTensor**](classTensor.md)&lt; VectorIndexSet&lt; [**R**](structR.md), [**Theta**](structTheta.md) &gt;, VectorIndexSet&lt; [**X\_cov**](classCircularToCartesian.md#typedef-x_cov), [**Y\_cov**](classCircularToCartesian.md#typedef-y_cov) &gt; &gt; | [**inv\_jacobian\_matrix**](#function-inv_jacobian_matrix) (Coord&lt; [**R**](structR.md), [**Theta**](structTheta.md) &gt; const & coord) const<br>_Compute full inverse Jacobian matrix._  |
|  KOKKOS\_FUNCTION double | [**jacobian**](#function-jacobian-12) (Coord&lt; [**R**](structR.md), [**Theta**](structTheta.md) &gt; const & coord) const<br>_Compute the Jacobian, the determinant of the Jacobian matrix of the mapping._  |
|  KOKKOS\_INLINE\_FUNCTION double | [**jacobian**](#function-jacobian-22) (Coord&lt; [**R**](structR.md) &gt; const & coord) const<br>_Compute the Jacobian, the determinant of the Jacobian matrix of the mapping as a function of the radius._  |
|  KOKKOS\_FUNCTION double | [**jacobian\_component**](#function-jacobian_component) (Coord&lt; [**R**](structR.md), [**Theta**](structTheta.md) &gt; const & coord) const<br>_Compute the (i,j) coefficient of the Jacobian matrix._  |
|  KOKKOS\_FUNCTION [**DTensor**](classTensor.md)&lt; VectorIndexSet&lt; [**X**](structX.md), [**Y**](structY.md) &gt;, VectorIndexSet&lt; [**R\_cov**](classCircularToCartesian.md#typedef-r_cov), [**Theta\_cov**](classCircularToCartesian.md#typedef-theta_cov) &gt; &gt; | [**jacobian\_matrix**](#function-jacobian_matrix) (Coord&lt; [**R**](structR.md), [**Theta**](structTheta.md) &gt; const & coord) const<br>_Compute full Jacobian matrix._  |
|  KOKKOS\_INLINE\_FUNCTION Coord&lt; [**X**](structX.md), [**Y**](structY.md) &gt; | [**o\_point**](#function-o_point) () const<br>_Get the O-point in Cartesian coordinates._  |
|  KOKKOS\_FUNCTION Coord&lt; [**X**](structX.md), [**Y**](structY.md) &gt; | [**operator()**](#function-operator) (Coord&lt; [**R**](structR.md), [**Theta**](structTheta.md) &gt; const & coord) const<br>_Convert the_ \((r, \theta)\) _coordinate to the equivalent (x,y) coordinate._ |
|  [**CircularToCartesian**](classCircularToCartesian.md) & | [**operator=**](#function-operator_1) ([**CircularToCartesian**](classCircularToCartesian.md) const & x) = default<br>_Assign a_ [_**CircularToCartesian**_](classCircularToCartesian.md) _from another_[_**CircularToCartesian**_](classCircularToCartesian.md) _(lvalue)._ |
|  [**CircularToCartesian**](classCircularToCartesian.md) & | [**operator=**](#function-operator_2) ([**CircularToCartesian**](classCircularToCartesian.md) && x) = default<br>_Assign a_ [_**CircularToCartesian**_](classCircularToCartesian.md) _from another temporary_[_**CircularToCartesian**_](classCircularToCartesian.md) _(rvalue)._ |
|  KOKKOS\_DEFAULTED\_FUNCTION | [**~CircularToCartesian**](#function-circulartocartesian) () = default<br> |




























## Detailed Description


The mapping \((r,\theta)\mapsto (x,y)\) is defined as follow :



\[x(r,\theta) = r \cos(\theta) + x_0,\]




\[y(r,\theta) = r \sin(\theta) + y_0.\]



It and its Jacobian matrix are invertible everywhere except for \(r = 0\).


The Jacobian matrix coefficients are defined as follow



\[J_{11}(r,\theta)  = \cos(\theta)\]




\[J_{12}(r,\theta)  = - r \sin(\theta)\]




\[J_{21}(r,\theta)  = \sin(\theta)\]




\[J_{22}(r,\theta)  = r \cos(\theta)\]



and the matrix determinant: \(det(J) = r\). 


    
## Public Types Documentation




### typedef CoordArg 

_The type of the argument of the function described by this mapping._ 
```C++
using CircularToCartesian< R, Theta, X, Y >::CoordArg =  Coord<R, Theta>;
```




<hr>



### typedef CoordJacobian 

_The type of the coordinate that can be used to evaluate the Jacobian of this mapping._ 
```C++
using CircularToCartesian< R, Theta, X, Y >::CoordJacobian =  CoordArg;
```




<hr>



### typedef CoordResult 

_The type of the result of the function described by this mapping._ 
```C++
using CircularToCartesian< R, Theta, X, Y >::CoordResult =  Coord<X, Y>;
```




<hr>



### typedef R\_cov 

_The covariant form of the first logical coordinate._ 
```C++
using CircularToCartesian< R, Theta, X, Y >::R_cov =  typename R::Dual;
```




<hr>



### typedef Theta\_cov 

_The covariant form of the second logical coordinate._ 
```C++
using CircularToCartesian< R, Theta, X, Y >::Theta_cov =  typename Theta::Dual;
```




<hr>



### typedef X\_cov 

_The covariant form of the first physical coordinate._ 
```C++
using CircularToCartesian< R, Theta, X, Y >::X_cov =  typename X::Dual;
```




<hr>



### typedef Y\_cov 

_The covariant form of the second physical coordinate._ 
```C++
using CircularToCartesian< R, Theta, X, Y >::Y_cov =  typename Y::Dual;
```




<hr>



### typedef cartesian\_tag\_x 

_Indicate the first physical coordinate._ 
```C++
using CircularToCartesian< R, Theta, X, Y >::cartesian_tag_x =  X;
```




<hr>



### typedef cartesian\_tag\_y 

_Indicate the second physical coordinate._ 
```C++
using CircularToCartesian< R, Theta, X, Y >::cartesian_tag_y =  Y;
```




<hr>



### typedef curvilinear\_tag\_r 

_Indicate the first logical coordinate._ 
```C++
using CircularToCartesian< R, Theta, X, Y >::curvilinear_tag_r =  R;
```




<hr>



### typedef curvilinear\_tag\_theta 

_Indicate the second logical coordinate._ 
```C++
using CircularToCartesian< R, Theta, X, Y >::curvilinear_tag_theta =  Theta;
```




<hr>
## Public Functions Documentation




### function CircularToCartesian [1/3]

_Instantiate a_ [_**CircularToCartesian**_](classCircularToCartesian.md) _from parameters._
```C++
inline explicit KOKKOS_FUNCTION CircularToCartesian::CircularToCartesian (
    Coord< X , Y > o_point=Coord< X , Y >(0.0, 0.0)
) 
```





**Parameters:**


* `o_point` The (x,y)-coordinate of the centre of the circle ((0,0) by default). 




        

<hr>



### function CircularToCartesian [2/3]

_Instantiate a_ [_**CircularToCartesian**_](classCircularToCartesian.md) _from another_[_**CircularToCartesian**_](classCircularToCartesian.md) _(lvalue)._
```C++
KOKKOS_DEFAULTED_FUNCTION CircularToCartesian::CircularToCartesian (
    CircularToCartesian const & other
) = default
```





**Parameters:**


* `other` [**CircularToCartesian**](classCircularToCartesian.md) mapping used to instantiate the new one. 




        

<hr>



### function CircularToCartesian [3/3]

_Instantiate a_ [_**CircularToCartesian**_](classCircularToCartesian.md) _from another temporary_[_**CircularToCartesian**_](classCircularToCartesian.md) _(rvalue)._
```C++
CircularToCartesian::CircularToCartesian (
    CircularToCartesian && x
) = default
```





**Parameters:**


* `x` [**CircularToCartesian**](classCircularToCartesian.md) mapping used to instantiate the new one. 




        

<hr>



### function get\_inverse\_mapping 

_Get the inverse mapping._ 
```C++
inline KOKKOS_INLINE_FUNCTION CartesianToCircular < X , Y , R , Theta > CircularToCartesian::get_inverse_mapping () const
```





**Returns:**

The inverse mapping. 





        

<hr>



### function inv\_jacobian\_component 

_Compute the (i,j) coefficient of the inverse Jacobian matrix._ 
```C++
template<class IndexTag1, class IndexTag2>
inline KOKKOS_FUNCTION double CircularToCartesian::inv_jacobian_component (
    Coord< R , Theta > const & coord
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
inline KOKKOS_FUNCTION DTensor < VectorIndexSet< R , Theta >, VectorIndexSet< X_cov , Y_cov > > CircularToCartesian::inv_jacobian_matrix (
    Coord< R , Theta > const & coord
) const
```



For some computations, we need the complete inverse Jacobian matrix or just the coefficients. The coefficients can be given independently with the function inv\_jacobian\_component.




**Parameters:**


* `coord` The coordinate where we evaluate the Jacobian matrix. 



**Returns:**

The inverse Jacobian matrix. 





        

<hr>



### function jacobian [1/2]

_Compute the Jacobian, the determinant of the Jacobian matrix of the mapping._ 
```C++
inline KOKKOS_FUNCTION double CircularToCartesian::jacobian (
    Coord< R , Theta > const & coord
) const
```





**Parameters:**


* `coord` The coordinate where we evaluate the Jacobian.



**Returns:**

A double with the value of the determinant of the Jacobian matrix. 





        

<hr>



### function jacobian [2/2]

_Compute the Jacobian, the determinant of the Jacobian matrix of the mapping as a function of the radius._ 
```C++
inline KOKKOS_INLINE_FUNCTION double CircularToCartesian::jacobian (
    Coord< R > const & coord
) const
```





**Parameters:**


* `coord` The radial coordinate where we evaluate the Jacobian.



**Returns:**

A double with the value of the determinant of the Jacobian matrix. 





        

<hr>



### function jacobian\_component 

_Compute the (i,j) coefficient of the Jacobian matrix._ 
```C++
template<class IndexTag1, class IndexTag2>
inline KOKKOS_FUNCTION double CircularToCartesian::jacobian_component (
    Coord< R , Theta > const & coord
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
inline KOKKOS_FUNCTION DTensor < VectorIndexSet< X , Y >, VectorIndexSet< R_cov , Theta_cov > > CircularToCartesian::jacobian_matrix (
    Coord< R , Theta > const & coord
) const
```



For some computations, we need the complete Jacobian matrix or just the coefficients. The coefficients can be given independently with the function jacobian\_component.




**Parameters:**


* `coord` The coordinate where we evaluate the Jacobian matrix. 



**Returns:**

The Jacobian matrix. 





        

<hr>



### function o\_point 

_Get the O-point in Cartesian coordinates._ 
```C++
inline KOKKOS_INLINE_FUNCTION Coord< X , Y > CircularToCartesian::o_point () const
```



Get the O-point of this mapping in Cartesian coordinates. This is calculated in the constructor.




**Returns:**

The O-point. 





        

<hr>



### function operator() 

_Convert the_ \((r, \theta)\) _coordinate to the equivalent (x,y) coordinate._
```C++
inline KOKKOS_FUNCTION Coord< X , Y > CircularToCartesian::operator() (
    Coord< R , Theta > const & coord
) const
```





**Parameters:**


* `coord` The coordinate to be converted.



**Returns:**

The equivalent coordinate. 





        

<hr>



### function operator= 

_Assign a_ [_**CircularToCartesian**_](classCircularToCartesian.md) _from another_[_**CircularToCartesian**_](classCircularToCartesian.md) _(lvalue)._
```C++
CircularToCartesian & CircularToCartesian::operator= (
    CircularToCartesian const & x
) = default
```





**Parameters:**


* `x` [**CircularToCartesian**](classCircularToCartesian.md) mapping used to assign.



**Returns:**

The [**CircularToCartesian**](classCircularToCartesian.md) assigned. 





        

<hr>



### function operator= 

_Assign a_ [_**CircularToCartesian**_](classCircularToCartesian.md) _from another temporary_[_**CircularToCartesian**_](classCircularToCartesian.md) _(rvalue)._
```C++
CircularToCartesian & CircularToCartesian::operator= (
    CircularToCartesian && x
) = default
```





**Parameters:**


* `x` [**CircularToCartesian**](classCircularToCartesian.md) mapping used to assign.



**Returns:**

The [**CircularToCartesian**](classCircularToCartesian.md) assigned. 





        

<hr>



### function ~CircularToCartesian 

```C++
KOKKOS_DEFAULTED_FUNCTION CircularToCartesian::~CircularToCartesian () = default
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/coord_transformations/circular_to_cartesian.hpp`

