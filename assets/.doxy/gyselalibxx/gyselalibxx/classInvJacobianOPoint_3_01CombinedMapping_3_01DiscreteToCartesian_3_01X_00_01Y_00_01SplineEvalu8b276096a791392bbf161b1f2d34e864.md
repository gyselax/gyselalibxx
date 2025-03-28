

# Class InvJacobianOPoint&lt; CombinedMapping&lt; DiscreteToCartesian&lt; X, Y, SplineEvaluator, R, Theta, MemorySpace &gt;, CartesianToCircular&lt; Xpc, Ypc, R, Theta &gt; &gt;, Coord&lt; R, Theta &gt; &gt;

**template &lt;class [**X**](structX.md), class [**Y**](structY.md), class SplineEvaluator, class [**R**](structR.md), class [**Theta**](structTheta.md), class MemorySpace, class Xpc, class Ypc&gt;**



[**ClassList**](annotated.md) **>** [**InvJacobianOPoint&lt; CombinedMapping&lt; DiscreteToCartesian&lt; X, Y, SplineEvaluator, R, Theta, MemorySpace &gt;, CartesianToCircular&lt; Xpc, Ypc, R, Theta &gt; &gt;, Coord&lt; R, Theta &gt; &gt;**](classInvJacobianOPoint_3_01CombinedMapping_3_01DiscreteToCartesian_3_01X_00_01Y_00_01SplineEvalu8b276096a791392bbf161b1f2d34e864.md)



[More...](#detailed-description)

* `#include <inv_jacobian_o_point.hpp>`





































## Public Functions

| Type | Name |
| ---: | :--- |
|  KOKKOS\_FUNCTION | [**InvJacobianOPoint**](#function-invjacobianopoint) ([**Mapping**](classCombinedMapping.md) const & mapping) <br>_The constructor of_ [_**InvJacobianOPoint**_](classInvJacobianOPoint.md) _._ |
|  KOKKOS\_FUNCTION [**DTensor**](classTensor.md)&lt; VectorIndexSet&lt; Xpc, Ypc &gt;, VectorIndexSet&lt; X\_cov, Y\_cov &gt; &gt; | [**operator()**](#function-operator) () const<br>_Compute the full inverse Jacobian matrix from a coordinate system (x\_pc, y\_pc) to a coordinate system (x, y) at the central point._  |




























## Detailed Description


A specialisation of [**InvJacobianOPoint**](classInvJacobianOPoint.md) for a combined mapping  where  is a discrete mapping from logical to physical, and  is an inverse circular mapping from physical to logical. The combined mapping  therefore maps from a physical domain  to a physical domain . 


    
## Public Functions Documentation




### function InvJacobianOPoint 

_The constructor of_ [_**InvJacobianOPoint**_](classInvJacobianOPoint.md) _._
```C++
inline explicit KOKKOS_FUNCTION InvJacobianOPoint< CombinedMapping< DiscreteToCartesian< X, Y, SplineEvaluator, R, Theta, MemorySpace >, CartesianToCircular< Xpc, Ypc, R, Theta > >, Coord< R, Theta > >::InvJacobianOPoint (
    Mapping const & mapping
) 
```





**Parameters:**


* `mapping` The mapping for which the inverse of the Jacobian is calculated. 




        

<hr>



### function operator() 

_Compute the full inverse Jacobian matrix from a coordinate system (x\_pc, y\_pc) to a coordinate system (x, y) at the central point._ 
```C++
inline KOKKOS_FUNCTION DTensor < VectorIndexSet< Xpc, Ypc >, VectorIndexSet< X_cov, Y_cov > > InvJacobianOPoint< CombinedMapping< DiscreteToCartesian< X, Y, SplineEvaluator, R, Theta, MemorySpace >, CartesianToCircular< Xpc, Ypc, R, Theta > >, Coord< R, Theta > >::operator() () const
```



The discrete mappings can be difficult to inverse especially at the central point. In case of non analytical invertible mapping, we can work in another domain called pseudo-Cartesian domain. In this domain, it is easier to inverse the Jacobian matrix. The idea is detailed in Edoardo Zoni's article (_Solving hyperbolic-elliptic problems on singular mapped disk-like domains with the method of characteristics and spline finite elements_, [https://doi.org/10.1016/j.jcp.2019.108889](https://doi.org/10.1016/j.jcp.2019.108889))


The current mapping maps from the logical domain to the physical domain . The pseudo-Cartesian domain is built with  and using the circular mapping :
* ,
* .




The pseudo-Cartesian domain is obtained by the composition of both mappings: . This new mapping is invertible and its inverse at the central point is given by
* 
* 
* 
* 




So the pseudo-Cartesian Jacobian matrix at the central point, , is obtained by inversing this matrix.




**Returns:**

The matrix evaluated at the central point.




**See also:** BslAdvection 


**See also:** AdvectionDomain 



        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/mapping/inv_jacobian_o_point.hpp`

