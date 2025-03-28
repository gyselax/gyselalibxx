

# Class InvJacobianOPoint&lt; CombinedMapping&lt; CircularToCartesian&lt; R, Theta, X, Y &gt;, CartesianToCircular&lt; Xpc, Ypc, R, Theta &gt; &gt;, Coord&lt; R, Theta &gt; &gt;

**template &lt;class [**X**](structX.md), class [**Y**](structY.md), class [**R**](structR.md), class [**Theta**](structTheta.md), class Xpc, class Ypc&gt;**



[**ClassList**](annotated.md) **>** [**InvJacobianOPoint&lt; CombinedMapping&lt; CircularToCartesian&lt; R, Theta, X, Y &gt;, CartesianToCircular&lt; Xpc, Ypc, R, Theta &gt; &gt;, Coord&lt; R, Theta &gt; &gt;**](classInvJacobianOPoint_3_01CombinedMapping_3_01CircularToCartesian_3_01R_00_01Theta_00_01X_00_01be6b75c3c69e2165a260584a5fd55276.md)



[More...](#detailed-description)

* `#include <inv_jacobian_o_point.hpp>`





































## Public Functions

| Type | Name |
| ---: | :--- |
|  KOKKOS\_FUNCTION | [**InvJacobianOPoint**](#function-invjacobianopoint) ([**Mapping**](classCombinedMapping.md) const & mapping) <br>_The constructor of_ [_**InvJacobianOPoint**_](classInvJacobianOPoint.md) _._ |
|  KOKKOS\_INLINE\_FUNCTION [**DTensor**](classTensor.md)&lt; VectorIndexSet&lt; Xpc, Ypc &gt;, VectorIndexSet&lt; X\_cov, Y\_cov &gt; &gt; | [**operator()**](#function-operator) () const<br>_Compute the full inverse Jacobian matrix from a coordinate system (x\_pc, y\_pc) to a coordinate system (x, y)._  |




























## Detailed Description


A specialisation of [**InvJacobianOPoint**](classInvJacobianOPoint.md) for a combined mapping  where  is a circular mapping from logical to physical, and  is an inverse circular mapping from physical to logical. The combined mapping  therefore maps from a physical domain  to a physical domain  (this mapping is equivalent to the identity). 


    
## Public Functions Documentation




### function InvJacobianOPoint 

_The constructor of_ [_**InvJacobianOPoint**_](classInvJacobianOPoint.md) _._
```C++
inline explicit KOKKOS_FUNCTION InvJacobianOPoint< CombinedMapping< CircularToCartesian< R, Theta, X, Y >, CartesianToCircular< Xpc, Ypc, R, Theta > >, Coord< R, Theta > >::InvJacobianOPoint (
    Mapping const & mapping
) 
```





**Parameters:**


* `mapping` The mapping for which the inverse of the Jacobian is calculated. 




        

<hr>



### function operator() 

_Compute the full inverse Jacobian matrix from a coordinate system (x\_pc, y\_pc) to a coordinate system (x, y)._ 
```C++
inline KOKKOS_INLINE_FUNCTION DTensor < VectorIndexSet< Xpc, Ypc >, VectorIndexSet< X_cov, Y_cov > > InvJacobianOPoint< CombinedMapping< CircularToCartesian< R, Theta, X, Y >, CartesianToCircular< Xpc, Ypc, R, Theta > >, Coord< R, Theta > >::operator() () const
```



Here, as , the Jacobian matrix of  is the identity matrix. So, the pseudo-Cartesian Jacobian matrix for a circular mapping is given by :
* 
* 
* 
* 






**Returns:**

The matrix evaluated at the central point. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/mapping/inv_jacobian_o_point.hpp`

