

# Class InvJacobianOPoint&lt; CombinedMapping&lt; CzarnyToCartesian&lt; R, Theta, X, Y &gt;, CartesianToCircular&lt; Xpc, Ypc, R, Theta &gt;, Coord&lt; R, Theta &gt; &gt;, Coord&lt; R, Theta &gt; &gt;

**template &lt;class [**X**](structX.md), class [**Y**](structY.md), class [**R**](structR.md), class [**Theta**](structTheta.md), class Xpc, class Ypc&gt;**



[**ClassList**](annotated.md) **>** [**InvJacobianOPoint&lt; CombinedMapping&lt; CzarnyToCartesian&lt; R, Theta, X, Y &gt;, CartesianToCircular&lt; Xpc, Ypc, R, Theta &gt;, Coord&lt; R, Theta &gt; &gt;, Coord&lt; R, Theta &gt; &gt;**](classInvJacobianOPoint_3_01CombinedMapping_3_01CzarnyToCartesian_3_01R_00_01Theta_00_01X_00_01Y_f284f6a7d72ad542b1021d394c9404b9.md)



[More...](#detailed-description)

* `#include <inv_jacobian_o_point.hpp>`





































## Public Functions

| Type | Name |
| ---: | :--- |
|  KOKKOS\_FUNCTION | [**InvJacobianOPoint**](#function-invjacobianopoint) ([**Mapping**](classCombinedMapping.md) const & mapping) <br>_The constructor of_ [_**InvJacobianOPoint**_](classInvJacobianOPoint.md) _._ |
|  KOKKOS\_INLINE\_FUNCTION [**DTensor**](classTensor.md)&lt; VectorIndexSet&lt; Xpc, Ypc &gt;, VectorIndexSet&lt; X\_cov, Y\_cov &gt; &gt; | [**operator()**](#function-operator) () const<br>_Compute the full inverse Jacobian matrix from a coordinate system (x\_pc, y\_pc) to a coordinate system (x, y) at the central point._  |




























## Detailed Description


A specialisation of [**InvJacobianOPoint**](classInvJacobianOPoint.md) for a combined mapping \(\mathcal{F} \circ \mathcal{G}\) where \(\mathcal{F}\) is a Czarny mapping from logical to physical, and \(\mathcal{G}\) is an inverse circular mapping from physical to logical. The combined mapping \(\mathcal{F} \circ \mathcal{G}\) therefore maps from a physical domain \((X_{pc}, Y_{pc})\) to a physical domain \((X, Y)\). 


    
## Public Functions Documentation




### function InvJacobianOPoint 

_The constructor of_ [_**InvJacobianOPoint**_](classInvJacobianOPoint.md) _._
```C++
inline explicit KOKKOS_FUNCTION InvJacobianOPoint< CombinedMapping< CzarnyToCartesian< R, Theta, X, Y >, CartesianToCircular< Xpc, Ypc, R, Theta >, Coord< R, Theta > >, Coord< R, Theta > >::InvJacobianOPoint (
    Mapping const & mapping
) 
```





**Parameters:**


* `mapping` The mapping for which the inverse of the Jacobian is calculated. 




        

<hr>



### function operator() 

_Compute the full inverse Jacobian matrix from a coordinate system (x\_pc, y\_pc) to a coordinate system (x, y) at the central point._ 
```C++
inline KOKKOS_INLINE_FUNCTION DTensor < VectorIndexSet< Xpc, Ypc >, VectorIndexSet< X_cov, Y_cov > > InvJacobianOPoint< CombinedMapping< CzarnyToCartesian< R, Theta, X, Y >, CartesianToCircular< Xpc, Ypc, R, Theta >, Coord< R, Theta > >, Coord< R, Theta > >::operator() () const
```



The pseudo-Cartesian Jacobian matrix for a Czarny mapping is given by :
* 
  \[(J_{\mathcal{F}}J_{\mathcal{G}}^{-1})^{-1}_{11}(0, \theta) = - \sqrt{1 + \varepsilon^2},\]

* 
  \[(J_{\mathcal{F}}J_{\mathcal{G}}^{-1})^{-1}_{12}(0, \theta) = 0,\]

* 
  \[(J_{\mathcal{F}}J_{\mathcal{G}}^{-1})^{-1}_{21}(0, \theta) = 0,\]

* 
  \[(J_{\mathcal{F}}J_{\mathcal{G}}^{-1})^{-1}_{22}(0, \theta) = \frac{2 - \sqrt{1 + \varepsilon^2}}{e \xi}.\]







**Returns:**

The matrix evaluated at the central point. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/coord_transformations/inv_jacobian_o_point.hpp`

