

# Class InvJacobianOPoint&lt; CombinedMapping&lt; DiscreteToCartesian&lt; X, Y, SplineEvaluator, R, Theta, MemorySpace &gt;, CartesianToCircular&lt; Xpc, Ypc, R, Theta &gt;, Coord&lt; R, Theta &gt; &gt;, Coord&lt; R, Theta &gt; &gt;

**template &lt;class [**X**](structX.md), class [**Y**](structY.md), class SplineEvaluator, class [**R**](structR.md), class [**Theta**](structTheta.md), class MemorySpace, class Xpc, class Ypc&gt;**



[**ClassList**](annotated.md) **>** [**InvJacobianOPoint&lt; CombinedMapping&lt; DiscreteToCartesian&lt; X, Y, SplineEvaluator, R, Theta, MemorySpace &gt;, CartesianToCircular&lt; Xpc, Ypc, R, Theta &gt;, Coord&lt; R, Theta &gt; &gt;, Coord&lt; R, Theta &gt; &gt;**](classInvJacobianOPoint_3_01CombinedMapping_3_01DiscreteToCartesian_3_01X_00_01Y_00_01SplineEvalu87e172e6ebb8e90a8cd02328541a469b.md)



[More...](#detailed-description)

* `#include <inv_jacobian_o_point.hpp>`





































## Public Functions

| Type | Name |
| ---: | :--- |
|  KOKKOS\_FUNCTION | [**InvJacobianOPoint**](#function-invjacobianopoint) ([**Mapping**](classCombinedMapping.md) const & mapping) <br>_The constructor of_ [_**InvJacobianOPoint**_](classInvJacobianOPoint.md) _._ |
|  KOKKOS\_FUNCTION [**DTensor**](classTensor.md)&lt; VectorIndexSet&lt; Xpc, Ypc &gt;, VectorIndexSet&lt; X\_cov, Y\_cov &gt; &gt; | [**operator()**](#function-operator) () const<br>_Compute the full inverse Jacobian matrix from a coordinate system (x\_pc, y\_pc) to a coordinate system (x, y) at the central point._  |




























## Detailed Description


A specialisation of [**InvJacobianOPoint**](classInvJacobianOPoint.md) for a combined mapping \(\mathcal{F} \circ \mathcal{G}\) where \(\mathcal{F}\) is a discrete mapping from logical to physical, and \(\mathcal{G}\) is an inverse circular mapping from physical to logical. The combined mapping \(\mathcal{F} \circ \mathcal{G}\) therefore maps from a physical domain \((X_{pc}, Y_{pc})\) to a physical domain \((X, Y)\). 


    
## Public Functions Documentation




### function InvJacobianOPoint 

_The constructor of_ [_**InvJacobianOPoint**_](classInvJacobianOPoint.md) _._
```C++
inline explicit KOKKOS_FUNCTION InvJacobianOPoint< CombinedMapping< DiscreteToCartesian< X, Y, SplineEvaluator, R, Theta, MemorySpace >, CartesianToCircular< Xpc, Ypc, R, Theta >, Coord< R, Theta > >, Coord< R, Theta > >::InvJacobianOPoint (
    Mapping const & mapping
) 
```





**Parameters:**


* `mapping` The mapping for which the inverse of the Jacobian is calculated. 




        

<hr>



### function operator() 

_Compute the full inverse Jacobian matrix from a coordinate system (x\_pc, y\_pc) to a coordinate system (x, y) at the central point._ 
```C++
inline KOKKOS_FUNCTION DTensor < VectorIndexSet< Xpc, Ypc >, VectorIndexSet< X_cov, Y_cov > > InvJacobianOPoint< CombinedMapping< DiscreteToCartesian< X, Y, SplineEvaluator, R, Theta, MemorySpace >, CartesianToCircular< Xpc, Ypc, R, Theta >, Coord< R, Theta > >, Coord< R, Theta > >::operator() () const
```



The discrete mappings can be difficult to inverse especially at the central point. In case of non analytical invertible mapping, we can work in another domain called pseudo-Cartesian domain. In this domain, it is easier to inverse the Jacobian matrix. The idea is detailed in Edoardo Zoni's article (_Solving hyperbolic-elliptic problems on singular mapped disk-like domains with the method of characteristics and spline finite elements_, [https://doi.org/10.1016/j.jcp.2019.108889](https://doi.org/10.1016/j.jcp.2019.108889))


The current mapping maps from the logical domain to the physical domain \(\mathcal{F}: (r,\theta) \mapsto (x, y)\). The pseudo-Cartesian domain is built with \(\mathcal{F}\) and using the circular mapping \(\mathcal{G}\):
* \(\mathcal{G}_{11}(r, \theta) = \cos(\theta),
     \qquad\quad \mathcal{G}_{12}(r, \theta) = \sin(\theta)\),
* \(\mathcal{G}_{21}(r, \theta) = -\frac{1}{r}\sin(\theta),
     \qquad\quad \mathcal{G}_{22}(r, \theta) = \frac{1}{r}\cos(\theta)\).




The pseudo-Cartesian domain is obtained by the composition of both mappings: \((\mathcal{F} \circ \mathcal{G}^{-1})^{-1}\). This new mapping is invertible and its inverse at the central point is given by
* 
\[(J_{\mathcal{F}}J_{\mathcal{G}}^{-1})_{11}(0, \theta) = \partial_r x (0, \theta) \cos(\theta)
             - \partial_{r \theta} x (0, \theta) \sin(\theta),\]

* 
\[(J_{\mathcal{F}}J_{\mathcal{G}}^{-1})_{12}(0, \theta) = \partial_r x (0, \theta) \sin(\theta)
             + \partial_{r \theta} x (0, \theta) \cos(\theta),\]

* 
\[(J_{\mathcal{F}}J_{\mathcal{G}}^{-1})_{21}(0, \theta) = \partial_r y (0, \theta) \cos(\theta)
             - \partial_{r \theta} y (0, \theta) \sin(\theta),\]

* 
\[(J_{\mathcal{F}}J_{\mathcal{G}}^{-1})_{22}(0, \theta) = \partial_r y (0, \theta) \sin(\theta)
             + \partial_{r \theta} y (0, \theta) \cos(\theta).\]





So the pseudo-Cartesian Jacobian matrix at the central point, \((J_{\mathcal{F}}J_{\mathcal{G}}^{-1})^{-1}(0, \theta)\), is obtained by inversing this matrix.




**Returns:**

The matrix evaluated at the central point.




**See also:** BslAdvection 


**See also:** AdvectionDomain 



        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/coord_transformations/inv_jacobian_o_point.hpp`

