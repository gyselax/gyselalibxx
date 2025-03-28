

# Class CzarnyToCartesian

**template &lt;class [**R**](structR.md), class [**Theta**](structTheta.md), class [**X**](structX.md), class [**Y**](structY.md)&gt;**



[**ClassList**](annotated.md) **>** [**CzarnyToCartesian**](classCzarnyToCartesian.md)



_A class for describing the Czarny 2D mapping._ [More...](#detailed-description)

* `#include <czarny_to_cartesian.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef Coord&lt; [**R**](structR.md), [**Theta**](structTheta.md) &gt; | [**CoordArg**](#typedef-coordarg)  <br>_The type of the argument of the function described by this mapping._  |
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
|   | [**CzarnyToCartesian**](#function-czarnytocartesian-13) (double epsilon, double e) <br>_Instantiate a_ [_**CzarnyToCartesian**_](classCzarnyToCartesian.md) _from parameters._ |
|  KOKKOS\_FUNCTION | [**CzarnyToCartesian**](#function-czarnytocartesian-23) ([**CzarnyToCartesian**](classCzarnyToCartesian.md) const & other) <br>_Instantiate a_ [_**CzarnyToCartesian**_](classCzarnyToCartesian.md) _from another_[_**CzarnyToCartesian**_](classCzarnyToCartesian.md) _(lvalue)._ |
|   | [**CzarnyToCartesian**](#function-czarnytocartesian-33) ([**CzarnyToCartesian**](classCzarnyToCartesian.md) && x) = default<br>_Instantiate a_ [_**CzarnyToCartesian**_](classCzarnyToCartesian.md) _from another temporary_[_**CzarnyToCartesian**_](classCzarnyToCartesian.md) _(rvalue)._ |
|  KOKKOS\_FUNCTION double | [**e**](#function-e) () const<br>_Return the_  _parameter._ |
|  KOKKOS\_FUNCTION double | [**epsilon**](#function-epsilon) () const<br>_Return the_  _parameter._ |
|  [**CartesianToCzarny**](classCartesianToCzarny.md)&lt; [**X**](structX.md), [**Y**](structY.md), [**R**](structR.md), [**Theta**](structTheta.md) &gt; | [**get\_inverse\_mapping**](#function-get_inverse_mapping) () const<br>_Get the inverse mapping._  |
|  KOKKOS\_FUNCTION double | [**inv\_jacobian\_11**](#function-inv_jacobian_11) (Coord&lt; [**R**](structR.md), [**Theta**](structTheta.md) &gt; const & coord) const<br>_Compute the (1,1) coefficient of the inverse Jacobian matrix._  |
|  KOKKOS\_FUNCTION double | [**inv\_jacobian\_12**](#function-inv_jacobian_12) (Coord&lt; [**R**](structR.md), [**Theta**](structTheta.md) &gt; const & coord) const<br>_Compute the (1,2) coefficient of the inverse Jacobian matrix._  |
|  KOKKOS\_FUNCTION double | [**inv\_jacobian\_21**](#function-inv_jacobian_21) (Coord&lt; [**R**](structR.md), [**Theta**](structTheta.md) &gt; const & coord) const<br>_Compute the (2,1) coefficient of the inverse Jacobian matrix._  |
|  KOKKOS\_FUNCTION double | [**inv\_jacobian\_22**](#function-inv_jacobian_22) (Coord&lt; [**R**](structR.md), [**Theta**](structTheta.md) &gt; const & coord) const<br>_Compute the (2,2) coefficient of the inverse Jacobian matrix._  |
|  KOKKOS\_FUNCTION [**DTensor**](classTensor.md)&lt; VectorIndexSet&lt; [**R**](structR.md), [**Theta**](structTheta.md) &gt;, VectorIndexSet&lt; [**X\_cov**](classCzarnyToCartesian.md#typedef-x_cov), [**Y\_cov**](classCzarnyToCartesian.md#typedef-y_cov) &gt; &gt; | [**inv\_jacobian\_matrix**](#function-inv_jacobian_matrix) (Coord&lt; [**R**](structR.md), [**Theta**](structTheta.md) &gt; const & coord) const<br>_Compute full inverse Jacobian matrix._  |
|  KOKKOS\_FUNCTION double | [**jacobian**](#function-jacobian) (Coord&lt; [**R**](structR.md), [**Theta**](structTheta.md) &gt; const & coord) const<br>_Compute the Jacobian, the determinant of the Jacobian matrix of the mapping._  |
|  KOKKOS\_FUNCTION double | [**jacobian\_component**](#function-jacobian_component) (Coord&lt; [**R**](structR.md), [**Theta**](structTheta.md) &gt; const & coord) const<br>_Compute the (i,j) coefficient of the Jacobian matrix._  |
|  KOKKOS\_FUNCTION [**DTensor**](classTensor.md)&lt; VectorIndexSet&lt; [**X**](structX.md), [**Y**](structY.md) &gt;, VectorIndexSet&lt; [**R\_cov**](classCzarnyToCartesian.md#typedef-r_cov), [**Theta\_cov**](classCzarnyToCartesian.md#typedef-theta_cov) &gt; &gt; | [**jacobian\_matrix**](#function-jacobian_matrix) (Coord&lt; [**R**](structR.md), [**Theta**](structTheta.md) &gt; const & coord) const<br>_Compute full Jacobian matrix._  |
|  KOKKOS\_FUNCTION Coord&lt; [**X**](structX.md), [**Y**](structY.md) &gt; | [**operator()**](#function-operator) (Coord&lt; [**R**](structR.md), [**Theta**](structTheta.md) &gt; const & coord) const<br>_Convert the_  _coordinate to the equivalent (x,y) coordinate._ |
|  [**CzarnyToCartesian**](classCzarnyToCartesian.md) & | [**operator=**](#function-operator_1) ([**CzarnyToCartesian**](classCzarnyToCartesian.md) const & x) = default<br>_Assign a_ [_**CzarnyToCartesian**_](classCzarnyToCartesian.md) _from another_[_**CzarnyToCartesian**_](classCzarnyToCartesian.md) _(lvalue)._ |
|  [**CzarnyToCartesian**](classCzarnyToCartesian.md) & | [**operator=**](#function-operator_2) ([**CzarnyToCartesian**](classCzarnyToCartesian.md) && x) = default<br>_Assign a_ [_**CzarnyToCartesian**_](classCzarnyToCartesian.md) _from another temporary_[_**CzarnyToCartesian**_](classCzarnyToCartesian.md) _(rvalue)._ |
|   | [**~CzarnyToCartesian**](#function-czarnytocartesian) () = default<br> |




























## Detailed Description


The mapping  is defined by








with  and  and  given as parameters. It and its Jacobian matrix are invertible everywhere except for .


Its Jacobian coefficients are defined as follow











.


and   


    
## Public Types Documentation




### typedef CoordArg 

_The type of the argument of the function described by this mapping._ 
```C++
using CzarnyToCartesian< R, Theta, X, Y >::CoordArg =  Coord<R, Theta>;
```




<hr>



### typedef CoordResult 

_The type of the result of the function described by this mapping._ 
```C++
using CzarnyToCartesian< R, Theta, X, Y >::CoordResult =  Coord<X, Y>;
```




<hr>



### typedef R\_cov 

_The covariant form of the first logical coordinate._ 
```C++
using CzarnyToCartesian< R, Theta, X, Y >::R_cov =  typename R::Dual;
```




<hr>



### typedef Theta\_cov 

_The covariant form of the second logical coordinate._ 
```C++
using CzarnyToCartesian< R, Theta, X, Y >::Theta_cov =  typename Theta::Dual;
```




<hr>



### typedef X\_cov 

_The covariant form of the first physical coordinate._ 
```C++
using CzarnyToCartesian< R, Theta, X, Y >::X_cov =  typename X::Dual;
```




<hr>



### typedef Y\_cov 

_The covariant form of the second physical coordinate._ 
```C++
using CzarnyToCartesian< R, Theta, X, Y >::Y_cov =  typename Y::Dual;
```




<hr>



### typedef cartesian\_tag\_x 

_Indicate the first physical coordinate._ 
```C++
using CzarnyToCartesian< R, Theta, X, Y >::cartesian_tag_x =  X;
```




<hr>



### typedef cartesian\_tag\_y 

_Indicate the second physical coordinate._ 
```C++
using CzarnyToCartesian< R, Theta, X, Y >::cartesian_tag_y =  Y;
```




<hr>



### typedef curvilinear\_tag\_r 

_Indicate the first logical coordinate._ 
```C++
using CzarnyToCartesian< R, Theta, X, Y >::curvilinear_tag_r =  R;
```




<hr>



### typedef curvilinear\_tag\_theta 

_Indicate the second logical coordinate._ 
```C++
using CzarnyToCartesian< R, Theta, X, Y >::curvilinear_tag_theta =  Theta;
```




<hr>
## Public Functions Documentation




### function CzarnyToCartesian [1/3]

_Instantiate a_ [_**CzarnyToCartesian**_](classCzarnyToCartesian.md) _from parameters._
```C++
inline CzarnyToCartesian::CzarnyToCartesian (
    double epsilon,
    double e
) 
```





**Parameters:**


* `epsilon` The  parameter in the definition of the mapping [**CzarnyToCartesian**](classCzarnyToCartesian.md).
* `e` The  parameter in the definition of the mapping [**CzarnyToCartesian**](classCzarnyToCartesian.md).



**See also:** [**CzarnyToCartesian**](classCzarnyToCartesian.md) 



        

<hr>



### function CzarnyToCartesian [2/3]

_Instantiate a_ [_**CzarnyToCartesian**_](classCzarnyToCartesian.md) _from another_[_**CzarnyToCartesian**_](classCzarnyToCartesian.md) _(lvalue)._
```C++
inline KOKKOS_FUNCTION CzarnyToCartesian::CzarnyToCartesian (
    CzarnyToCartesian const & other
) 
```





**Parameters:**


* `other` [**CzarnyToCartesian**](classCzarnyToCartesian.md) mapping used to instantiate the new one. 




        

<hr>



### function CzarnyToCartesian [3/3]

_Instantiate a_ [_**CzarnyToCartesian**_](classCzarnyToCartesian.md) _from another temporary_[_**CzarnyToCartesian**_](classCzarnyToCartesian.md) _(rvalue)._
```C++
CzarnyToCartesian::CzarnyToCartesian (
    CzarnyToCartesian && x
) = default
```





**Parameters:**


* `x` [**CzarnyToCartesian**](classCzarnyToCartesian.md) mapping used to instantiate the new one. 




        

<hr>



### function e 

_Return the_  _parameter._
```C++
inline KOKKOS_FUNCTION double CzarnyToCartesian::e () const
```





**Returns:**

The value of .




**See also:** [**CzarnyToCartesian**](classCzarnyToCartesian.md) 



        

<hr>



### function epsilon 

_Return the_  _parameter._
```C++
inline KOKKOS_FUNCTION double CzarnyToCartesian::epsilon () const
```





**Returns:**

The value of .




**See also:** [**CzarnyToCartesian**](classCzarnyToCartesian.md) 



        

<hr>



### function get\_inverse\_mapping 

_Get the inverse mapping._ 
```C++
inline CartesianToCzarny < X , Y , R , Theta > CzarnyToCartesian::get_inverse_mapping () const
```





**Returns:**

The inverse mapping. 





        

<hr>



### function inv\_jacobian\_11 

_Compute the (1,1) coefficient of the inverse Jacobian matrix._ 
```C++
inline KOKKOS_FUNCTION double CzarnyToCartesian::inv_jacobian_11 (
    Coord< R , Theta > const & coord
) const
```



Be careful because not all mappings are invertible, especially at the centre point.




**Parameters:**


* `coord` The coordinate where we evaluate the inverse Jacobian matrix.



**Returns:**

A double with the value of the (1,1) coefficient of the inverse Jacobian matrix. 





        

<hr>



### function inv\_jacobian\_12 

_Compute the (1,2) coefficient of the inverse Jacobian matrix._ 
```C++
inline KOKKOS_FUNCTION double CzarnyToCartesian::inv_jacobian_12 (
    Coord< R , Theta > const & coord
) const
```



Be careful because not all mappings are invertible, especially at the centre point.




**Parameters:**


* `coord` The coordinate where we evaluate the inverse Jacobian matrix.



**Returns:**

A double with the value of the (1,2) coefficient of the inverse Jacobian matrix. 





        

<hr>



### function inv\_jacobian\_21 

_Compute the (2,1) coefficient of the inverse Jacobian matrix._ 
```C++
inline KOKKOS_FUNCTION double CzarnyToCartesian::inv_jacobian_21 (
    Coord< R , Theta > const & coord
) const
```



Be careful because not all mappings are invertible, especially at the centre point.




**Parameters:**


* `coord` The coordinate where we evaluate the inverse Jacobian matrix.



**Returns:**

A double with the value of the (2,1) coefficient of the inverse Jacobian matrix. 





        

<hr>



### function inv\_jacobian\_22 

_Compute the (2,2) coefficient of the inverse Jacobian matrix._ 
```C++
inline KOKKOS_FUNCTION double CzarnyToCartesian::inv_jacobian_22 (
    Coord< R , Theta > const & coord
) const
```



Be careful because not all mappings are invertible, especially at the centre point.




**Parameters:**


* `coord` The coordinate where we evaluate the inverse Jacobian matrix.



**Returns:**

A double with the value of the (2,2) coefficient of the inverse Jacobian matrix. 





        

<hr>



### function inv\_jacobian\_matrix 

_Compute full inverse Jacobian matrix._ 
```C++
inline KOKKOS_FUNCTION DTensor < VectorIndexSet< R , Theta >, VectorIndexSet< X_cov , Y_cov > > CzarnyToCartesian::inv_jacobian_matrix (
    Coord< R , Theta > const & coord
) const
```



For some computations, we need the complete inverse Jacobian matrix or just the coefficients. The coefficients can be given independently with the functions inv\_jacobian\_11, inv\_jacobian\_12, inv\_jacobian\_21 and inv\_jacobian\_22.




**Parameters:**


* `coord` The coordinate where we evaluate the Jacobian matrix. 



**Returns:**

The inverse Jacobian matrix.




**See also:** Jacobian::inv\_jacobian\_11 


**See also:** Jacobian::inv\_jacobian\_12 


**See also:** Jacobian::inv\_jacobian\_21 


**See also:** Jacobian::inv\_jacobian\_22 



        

<hr>



### function jacobian 

_Compute the Jacobian, the determinant of the Jacobian matrix of the mapping._ 
```C++
inline KOKKOS_FUNCTION double CzarnyToCartesian::jacobian (
    Coord< R , Theta > const & coord
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
inline KOKKOS_FUNCTION double CzarnyToCartesian::jacobian_component (
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
inline KOKKOS_FUNCTION DTensor < VectorIndexSet< X , Y >, VectorIndexSet< R_cov , Theta_cov > > CzarnyToCartesian::jacobian_matrix (
    Coord< R , Theta > const & coord
) const
```



For some computations, we need the complete Jacobian matrix or just the coefficients. The coefficients can be given independently with the function jacobian\_component




**Parameters:**


* `coord` The coordinate where we evaluate the Jacobian matrix. 



**Returns:**

The Jacobian matrix. 





        

<hr>



### function operator() 

_Convert the_  _coordinate to the equivalent (x,y) coordinate._
```C++
inline KOKKOS_FUNCTION Coord< X , Y > CzarnyToCartesian::operator() (
    Coord< R , Theta > const & coord
) const
```





**Parameters:**


* `coord` The coordinate to be converted.



**Returns:**

The equivalent coordinate. 





        

<hr>



### function operator= 

_Assign a_ [_**CzarnyToCartesian**_](classCzarnyToCartesian.md) _from another_[_**CzarnyToCartesian**_](classCzarnyToCartesian.md) _(lvalue)._
```C++
CzarnyToCartesian & CzarnyToCartesian::operator= (
    CzarnyToCartesian const & x
) = default
```





**Parameters:**


* `x` [**CzarnyToCartesian**](classCzarnyToCartesian.md) mapping used to assign.



**Returns:**

The [**CzarnyToCartesian**](classCzarnyToCartesian.md) assigned. 





        

<hr>



### function operator= 

_Assign a_ [_**CzarnyToCartesian**_](classCzarnyToCartesian.md) _from another temporary_[_**CzarnyToCartesian**_](classCzarnyToCartesian.md) _(rvalue)._
```C++
CzarnyToCartesian & CzarnyToCartesian::operator= (
    CzarnyToCartesian && x
) = default
```





**Parameters:**


* `x` [**CzarnyToCartesian**](classCzarnyToCartesian.md) mapping used to assign.



**Returns:**

The [**CzarnyToCartesian**](classCzarnyToCartesian.md) assigned. 





        

<hr>



### function ~CzarnyToCartesian 

```C++
CzarnyToCartesian::~CzarnyToCartesian () = default
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/mapping/czarny_to_cartesian.hpp`

