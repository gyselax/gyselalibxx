

# Class DiscreteToCartesian

**template &lt;class [**X**](structX.md), class [**Y**](structY.md), class SplineEvaluator, class [**R**](structR.md), class [**Theta**](structTheta.md), class MemorySpace&gt;**



[**ClassList**](annotated.md) **>** [**DiscreteToCartesian**](classDiscreteToCartesian.md)



_A class for describing discrete 2D mappings from the logical domain to the physical domain._ [More...](#detailed-description)

* `#include <discrete_to_cartesian.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef typename SplineEvaluator::bsplines\_type1 | [**BSplineR**](#typedef-bspliner)  <br>_Indicate the bspline type of the first logical dimension._  |
| typedef typename SplineEvaluator::bsplines\_type2 | [**BSplineTheta**](#typedef-bsplinetheta)  <br>_Indicate the bspline type of the second logical dimension._  |
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
|  KOKKOS\_FUNCTION | [**DiscreteToCartesian**](#function-discretetocartesian) (SplineType curvilinear\_to\_x, SplineType curvilinear\_to\_y, SplineEvaluator const & evaluator, IdxRangeRTheta idx\_range\_singular\_point) <br>_Instantiate a_ [_**DiscreteToCartesian**_](classDiscreteToCartesian.md) _from the coefficients of 2D splines approximating the mapping._ |
|  KOKKOS\_INLINE\_FUNCTION const Coord&lt; [**X**](structX.md), [**Y**](structY.md) &gt; | [**control\_point**](#function-control_point) (Idx&lt; [**BSplineR**](classDiscreteToCartesian.md#typedef-bspliner), [**BSplineTheta**](classDiscreteToCartesian.md#typedef-bsplinetheta) &gt; const & el) const<br>_Get a control point of the mapping on B-splines._  |
|  KOKKOS\_INLINE\_FUNCTION [**DTensor**](classTensor.md)&lt; VectorIndexSet&lt; [**X**](structX.md), [**Y**](structY.md) &gt;, VectorIndexSet&lt; [**R\_cov**](classDiscreteToCartesian.md#typedef-r_cov), [**Theta\_cov**](classDiscreteToCartesian.md#typedef-theta_cov) &gt; &gt; | [**first\_order\_jacobian\_matrix\_r\_rtheta**](#function-first_order_jacobian_matrix_r_rtheta) (Coord&lt; [**curvilinear\_tag\_r**](classDiscreteToCartesian.md#typedef-curvilinear_tag_r), [**curvilinear\_tag\_theta**](classDiscreteToCartesian.md#typedef-curvilinear_tag_theta) &gt; const & coord) const<br>_Get the first order expansion of the Jacobian matrix with the theta component divided by r. The expansion is carried out around_  _. The returned matrix_ _is defined as:_ __ __ __ _._ |
|  KOKKOS\_INLINE\_FUNCTION IdxRangeRTheta | [**idx\_range\_singular\_point**](#function-idx_range_singular_point) () const<br>_Get the index range describing the points which should be used to evaluate functions at the central point._  |
|  KOKKOS\_FUNCTION double | [**jacobian**](#function-jacobian) (Coord&lt; [**curvilinear\_tag\_r**](classDiscreteToCartesian.md#typedef-curvilinear_tag_r), [**curvilinear\_tag\_theta**](classDiscreteToCartesian.md#typedef-curvilinear_tag_theta) &gt; const & coord) const<br>_Compute the Jacobian, the determinant of the Jacobian matrix of the mapping._  |
|  KOKKOS\_INLINE\_FUNCTION double | [**jacobian\_component**](#function-jacobian_component) (Coord&lt; [**R**](structR.md), [**Theta**](structTheta.md) &gt; coord) const<br>_Compute the (i,j) coefficient of the Jacobian matrix._  |
|  KOKKOS\_FUNCTION [**DTensor**](classTensor.md)&lt; VectorIndexSet&lt; [**X**](structX.md), [**Y**](structY.md) &gt;, VectorIndexSet&lt; [**R\_cov**](classDiscreteToCartesian.md#typedef-r_cov), [**Theta\_cov**](classDiscreteToCartesian.md#typedef-theta_cov) &gt; &gt; | [**jacobian\_matrix**](#function-jacobian_matrix) (Coord&lt; [**R**](structR.md), [**Theta**](structTheta.md) &gt; const & coord) const<br>_Compute full Jacobian matrix._  |
|  KOKKOS\_FUNCTION Coord&lt; [**X**](structX.md), [**Y**](structY.md) &gt; | [**operator()**](#function-operator) (Coord&lt; [**curvilinear\_tag\_r**](classDiscreteToCartesian.md#typedef-curvilinear_tag_r), [**curvilinear\_tag\_theta**](classDiscreteToCartesian.md#typedef-curvilinear_tag_theta) &gt; const & coord) const<br>_Compute the physical coordinates from the logical coordinates._  |




























## Detailed Description


The mapping describe here is only defined on a grid. The [**DiscreteToCartesian**](classDiscreteToCartesian.md) class decomposes the mapping on B-splines to evaluate it on the physical domain.








This mapping could be costly to inverse. 


    
## Public Types Documentation




### typedef BSplineR 

_Indicate the bspline type of the first logical dimension._ 
```C++
using DiscreteToCartesian< X, Y, SplineEvaluator, R, Theta, MemorySpace >::BSplineR =  typename SplineEvaluator::bsplines_type1;
```




<hr>



### typedef BSplineTheta 

_Indicate the bspline type of the second logical dimension._ 
```C++
using DiscreteToCartesian< X, Y, SplineEvaluator, R, Theta, MemorySpace >::BSplineTheta =  typename SplineEvaluator::bsplines_type2;
```




<hr>



### typedef CoordArg 

_The type of the argument of the function described by this mapping._ 
```C++
using DiscreteToCartesian< X, Y, SplineEvaluator, R, Theta, MemorySpace >::CoordArg =  Coord<R, Theta>;
```




<hr>



### typedef CoordResult 

_The type of the result of the function described by this mapping._ 
```C++
using DiscreteToCartesian< X, Y, SplineEvaluator, R, Theta, MemorySpace >::CoordResult =  Coord<X, Y>;
```




<hr>



### typedef R\_cov 

_The covariant form of the first logical coordinate._ 
```C++
using DiscreteToCartesian< X, Y, SplineEvaluator, R, Theta, MemorySpace >::R_cov =  typename R::Dual;
```




<hr>



### typedef Theta\_cov 

_The covariant form of the second logical coordinate._ 
```C++
using DiscreteToCartesian< X, Y, SplineEvaluator, R, Theta, MemorySpace >::Theta_cov =  typename Theta::Dual;
```




<hr>



### typedef X\_cov 

_The covariant form of the first physical coordinate._ 
```C++
using DiscreteToCartesian< X, Y, SplineEvaluator, R, Theta, MemorySpace >::X_cov =  typename X::Dual;
```




<hr>



### typedef Y\_cov 

_The covariant form of the second physical coordinate._ 
```C++
using DiscreteToCartesian< X, Y, SplineEvaluator, R, Theta, MemorySpace >::Y_cov =  typename Y::Dual;
```




<hr>



### typedef cartesian\_tag\_x 

_Indicate the first physical coordinate._ 
```C++
using DiscreteToCartesian< X, Y, SplineEvaluator, R, Theta, MemorySpace >::cartesian_tag_x =  X;
```




<hr>



### typedef cartesian\_tag\_y 

_Indicate the second physical coordinate._ 
```C++
using DiscreteToCartesian< X, Y, SplineEvaluator, R, Theta, MemorySpace >::cartesian_tag_y =  Y;
```




<hr>



### typedef curvilinear\_tag\_r 

_Indicate the first logical coordinate._ 
```C++
using DiscreteToCartesian< X, Y, SplineEvaluator, R, Theta, MemorySpace >::curvilinear_tag_r =  R;
```




<hr>



### typedef curvilinear\_tag\_theta 

_Indicate the second logical coordinate._ 
```C++
using DiscreteToCartesian< X, Y, SplineEvaluator, R, Theta, MemorySpace >::curvilinear_tag_theta =  Theta;
```




<hr>
## Public Functions Documentation




### function DiscreteToCartesian 

_Instantiate a_ [_**DiscreteToCartesian**_](classDiscreteToCartesian.md) _from the coefficients of 2D splines approximating the mapping._
```C++
inline KOKKOS_FUNCTION DiscreteToCartesian::DiscreteToCartesian (
    SplineType curvilinear_to_x,
    SplineType curvilinear_to_y,
    SplineEvaluator const & evaluator,
    IdxRangeRTheta idx_range_singular_point
) 
```



A discrete mapping is a mapping whose values are known only at the mesh points of the grid. To interpolate the mapping, we use B-splines. The [**DiscreteToCartesian**](classDiscreteToCartesian.md) mapping is initialised from the coefficients in front of the basis splines which arise when we approximate the functions , and  (with  and  the physical dimensions in the logical domain) with Splines (using SplineBuilder2D). Then to interpolate the mapping, we will evaluate the decomposed functions on B-splines (see DiscreteToCartesian::operator()).


Here, the evaluator is given as input.




**Parameters:**


* `curvilinear_to_x` Bsplines coefficients of the first physical dimension in the logical domain.
* `curvilinear_to_y` Bsplines coefficients of the second physical dimension in the logical domain.
* `evaluator` The evaluator used to evaluate the mapping.
* `idx_range_singular_point` The index range describing the points which should be used to evaluate functions at the central point.



**See also:** SplineBuilder2D 


**See also:** DiscreteToCartesian::operator() 


**See also:** SplineBoundaryValue 



        

<hr>



### function control\_point 

_Get a control point of the mapping on B-splines._ 
```C++
inline KOKKOS_INLINE_FUNCTION const Coord< X , Y > DiscreteToCartesian::control_point (
    Idx< BSplineR , BSplineTheta > const & el
) const
```



The mapping  decomposed on B-splines can be identified by its control points  where  and  are the B-splines coefficients:


,


,


where  is the number of B-splines.


The control points can be obtained by interpolating the mapping on interpolation points (see GrevilleInterpolationPoints or KnotsAsInterpolationPoints). We can also note that the first control points  are equal to the pole ,  where  .




**Parameters:**


* `el` The number of the control point.



**Returns:**

The el-th control point.




**See also:** GrevilleInterpolationPoints 


**See also:** KnotsAsInterpolationPoints 



        

<hr>



### function first\_order\_jacobian\_matrix\_r\_rtheta 

_Get the first order expansion of the Jacobian matrix with the theta component divided by r. The expansion is carried out around_  _. The returned matrix_ _is defined as:_ __ __ __ _._
```C++
inline KOKKOS_INLINE_FUNCTION DTensor < VectorIndexSet< X , Y >, VectorIndexSet< R_cov , Theta_cov > > DiscreteToCartesian::first_order_jacobian_matrix_r_rtheta (
    Coord< curvilinear_tag_r , curvilinear_tag_theta > const & coord
) const
```





**Parameters:**


* `coord` The coordinate where we evaluate the Jacobian.



**Returns:**

The first order expansion of the Jacobian matrix with the theta component divided by r. 





        

<hr>



### function idx\_range\_singular\_point 

_Get the index range describing the points which should be used to evaluate functions at the central point._ 
```C++
inline KOKKOS_INLINE_FUNCTION IdxRangeRTheta DiscreteToCartesian::idx_range_singular_point () const
```





**Returns:**

An index range covering the O-point. 





        

<hr>



### function jacobian 

_Compute the Jacobian, the determinant of the Jacobian matrix of the mapping._ 
```C++
inline KOKKOS_FUNCTION double DiscreteToCartesian::jacobian (
    Coord< curvilinear_tag_r , curvilinear_tag_theta > const & coord
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
inline KOKKOS_INLINE_FUNCTION double DiscreteToCartesian::jacobian_component (
    Coord< R , Theta > coord
) const
```



For a mapping given by , with  the curvilinear coordinates and  the Cartesian coordinates, the (i,j) coefficient of the Jacobian matrix is given by .


As the mapping is decomposed on B-splines, it means it computes the derivatives of B-splines  (the derivatives are implemented in SplineEvaluator2D).




**Parameters:**


* `coord` The coordinate where we evaluate the Jacobian matrix.



**Returns:**

A double with the value of the (i,j) coefficient of the Jacobian matrix.




**See also:** SplineEvaluator2D 



        

<hr>



### function jacobian\_matrix 

_Compute full Jacobian matrix._ 
```C++
inline KOKKOS_FUNCTION DTensor < VectorIndexSet< X , Y >, VectorIndexSet< R_cov , Theta_cov > > DiscreteToCartesian::jacobian_matrix (
    Coord< R , Theta > const & coord
) const
```



For some computations, we need the complete Jacobian matrix or just the coefficients. The coefficients can be given independently with the function jacobian\_component.




**Parameters:**


* `coord` The coordinate where we evaluate the Jacobian matrix. 



**Returns:**

The Jacobian matrix. 





        

<hr>



### function operator() 

_Compute the physical coordinates from the logical coordinates._ 
```C++
inline KOKKOS_FUNCTION Coord< X , Y > DiscreteToCartesian::operator() (
    Coord< curvilinear_tag_r , curvilinear_tag_theta > const & coord
) const
```



It evaluates the decomposed mapping on B-splines at the coordinate point with a SplineEvaluator2D.




**Parameters:**


* `coord` The coordinates in the logical domain.



**Returns:**

The coordinates of the mapping in the physical domain.




**See also:** SplineEvaluator2D 



        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/mapping/discrete_to_cartesian.hpp`

