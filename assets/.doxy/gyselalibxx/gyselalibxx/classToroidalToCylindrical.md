

# Class ToroidalToCylindrical

**template &lt;class Curvilinear2DToCartesian, class Zeta, class Phi&gt;**



[**ClassList**](annotated.md) **>** [**ToroidalToCylindrical**](classToroidalToCylindrical.md)



_A class describing a coordinate change from a toroidal system of coordinates to a cylindrical system of coordinates. The toroidal coordinates are described by a polar plane_  _and a perpendicular dimension_ _. The cylindrical coordinates are_ _._ _describe a Cartesian slice._ _are therefore defined from this slice with a 2D coordinate change operator._ _is chosen to be equal to_ _to preserve the orientation of the axes (following the right-hand rule)._[More...](#detailed-description)

* `#include <toroidal_to_cylindrical.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef Coord&lt; Rho, [**Theta**](structTheta.md), Phi &gt; | [**CoordArg**](#typedef-coordarg)  <br>_The type of the argument of the function described by this mapping._  |
| typedef Coord&lt; [**R**](structR.md), Z, Zeta &gt; | [**CoordResult**](#typedef-coordresult)  <br>_The type of the result of the function described by this mapping._  |
| typedef typename Curvilinear2DToCartesian::cartesian\_tag\_x | [**cylindrical\_tag\_R**](#typedef-cylindrical_tag_r)  <br>_Indicate the first physical coordinate._  |
| typedef typename Curvilinear2DToCartesian::cartesian\_tag\_y | [**cylindrical\_tag\_Z**](#typedef-cylindrical_tag_z)  <br>_Indicate the second physical coordinate._  |
| typedef Zeta | [**cylindrical\_tag\_Zeta**](#typedef-cylindrical_tag_zeta)  <br>_Indicate the third physical coordinate._  |
| typedef Phi | [**toroidal\_tag\_phi**](#typedef-toroidal_tag_phi)  <br>_Indicate the third logical coordinate._  |
| typedef typename Curvilinear2DToCartesian::curvilinear\_tag\_r | [**toroidal\_tag\_rho**](#typedef-toroidal_tag_rho)  <br>_Indicate the first logical coordinate._  |
| typedef typename Curvilinear2DToCartesian::curvilinear\_tag\_theta | [**toroidal\_tag\_theta**](#typedef-toroidal_tag_theta)  <br>_Indicate the second logical coordinate._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**ToroidalToCylindrical**](#function-toroidaltocylindrical) (Curvilinear2DToCartesian const & mapping\_2d) <br>_Instantiate a_ [_**ToroidalToCylindrical**_](classToroidalToCylindrical.md) _mapping from a 2D Curvilinear2DToCartesian mapping._ |
|  Curvilinear2DToCartesian | [**get\_2d\_polar\_mapping**](#function-get_2d_polar_mapping) () const<br>_Get the mapping describing the polar plane._  |
|  KOKKOS\_INLINE\_FUNCTION double | [**inv\_jacobian\_component**](#function-inv_jacobian_component) ([**CoordArg**](classToroidalToCylindrical.md#typedef-coordarg) const & coord) const<br>_Compute the (i,j) coefficient of the inverse of the Jacobian matrix._  |
|  KOKKOS\_FUNCTION [**DTensor**](classTensor.md)&lt; VectorIndexSet&lt; Rho, [**Theta**](structTheta.md), Phi &gt;, VectorIndexSet&lt; [**R\_cov**](structR__cov.md), Z\_cov, Zeta\_cov &gt; &gt; | [**inv\_jacobian\_matrix**](#function-inv_jacobian_matrix) ([**CoordArg**](classToroidalToCylindrical.md#typedef-coordarg) const & coord) const<br>_Compute the inverse of the Jacobian matrix._  |
|  KOKKOS\_FUNCTION double | [**jacobian**](#function-jacobian) ([**CoordArg**](classToroidalToCylindrical.md#typedef-coordarg) const & coord) const<br>_Compute the Jacobian: the determinant of the Jacobian matrix of the mapping._  |
|  KOKKOS\_INLINE\_FUNCTION double | [**jacobian\_component**](#function-jacobian_component) ([**CoordArg**](classToroidalToCylindrical.md#typedef-coordarg) const & coord) const<br>_Compute the (i,j) coefficient of the Jacobian matrix._  |
|  KOKKOS\_FUNCTION [**DTensor**](classTensor.md)&lt; VectorIndexSet&lt; [**R**](structR.md), Z, Zeta &gt;, VectorIndexSet&lt; Rho\_cov, [**Theta\_cov**](structTheta__cov.md), Phi\_cov &gt; &gt; | [**jacobian\_matrix**](#function-jacobian_matrix) ([**CoordArg**](classToroidalToCylindrical.md#typedef-coordarg) const & coord) const<br>_Compute the Jacobian matrix._  |
|  KOKKOS\_FUNCTION [**CoordResult**](classToroidalToCylindrical.md#typedef-coordresult) | [**operator()**](#function-operator) ([**CoordArg**](classToroidalToCylindrical.md#typedef-coordarg) const & coord) const<br>_Convert the_  _coordinate to the equivalent_ _coordinate._ |




























## Detailed Description




**Template parameters:**


* `Curvilinear2DToCartesian` An operator describing the coordinate change from  to . 
* `Zeta` The angle of the cylindrical coordinates. 
* `Phi` The toroidal component of the toroidal coordinates. 




    
## Public Types Documentation




### typedef CoordArg 

_The type of the argument of the function described by this mapping._ 
```C++
using ToroidalToCylindrical< Curvilinear2DToCartesian, Zeta, Phi >::CoordArg =  Coord<Rho, Theta, Phi>;
```




<hr>



### typedef CoordResult 

_The type of the result of the function described by this mapping._ 
```C++
using ToroidalToCylindrical< Curvilinear2DToCartesian, Zeta, Phi >::CoordResult =  Coord<R, Z, Zeta>;
```




<hr>



### typedef cylindrical\_tag\_R 

_Indicate the first physical coordinate._ 
```C++
using ToroidalToCylindrical< Curvilinear2DToCartesian, Zeta, Phi >::cylindrical_tag_R =  typename Curvilinear2DToCartesian::cartesian_tag_x;
```




<hr>



### typedef cylindrical\_tag\_Z 

_Indicate the second physical coordinate._ 
```C++
using ToroidalToCylindrical< Curvilinear2DToCartesian, Zeta, Phi >::cylindrical_tag_Z =  typename Curvilinear2DToCartesian::cartesian_tag_y;
```




<hr>



### typedef cylindrical\_tag\_Zeta 

_Indicate the third physical coordinate._ 
```C++
using ToroidalToCylindrical< Curvilinear2DToCartesian, Zeta, Phi >::cylindrical_tag_Zeta =  Zeta;
```




<hr>



### typedef toroidal\_tag\_phi 

_Indicate the third logical coordinate._ 
```C++
using ToroidalToCylindrical< Curvilinear2DToCartesian, Zeta, Phi >::toroidal_tag_phi =  Phi;
```




<hr>



### typedef toroidal\_tag\_rho 

_Indicate the first logical coordinate._ 
```C++
using ToroidalToCylindrical< Curvilinear2DToCartesian, Zeta, Phi >::toroidal_tag_rho =  typename Curvilinear2DToCartesian::curvilinear_tag_r;
```




<hr>



### typedef toroidal\_tag\_theta 

_Indicate the second logical coordinate._ 
```C++
using ToroidalToCylindrical< Curvilinear2DToCartesian, Zeta, Phi >::toroidal_tag_theta =  typename Curvilinear2DToCartesian::curvilinear_tag_theta;
```




<hr>
## Public Functions Documentation




### function ToroidalToCylindrical 

_Instantiate a_ [_**ToroidalToCylindrical**_](classToroidalToCylindrical.md) _mapping from a 2D Curvilinear2DToCartesian mapping._
```C++
inline explicit ToroidalToCylindrical::ToroidalToCylindrical (
    Curvilinear2DToCartesian const & mapping_2d
) 
```





**Parameters:**


* `mapping_2d` The mapping governing the transformation from  to . 




        

<hr>



### function get\_2d\_polar\_mapping 

_Get the mapping describing the polar plane._ 
```C++
inline Curvilinear2DToCartesian ToroidalToCylindrical::get_2d_polar_mapping () const
```





**Returns:**

The mapping describing the polar plane. 





        

<hr>



### function inv\_jacobian\_component 

_Compute the (i,j) coefficient of the inverse of the Jacobian matrix._ 
```C++
template<class IndexTag1, class IndexTag2>
inline KOKKOS_INLINE_FUNCTION double ToroidalToCylindrical::inv_jacobian_component (
    CoordArg const & coord
) const
```





**Parameters:**


* `coord` The coordinate where we evaluate the Jacobian matrix.



**Template parameters:**


* `The` tag representing the index i. 
* `The` tag representing the index j.



**Returns:**

The value of the (i,j) coefficient of the Jacobian matrix. 





        

<hr>



### function inv\_jacobian\_matrix 

_Compute the inverse of the Jacobian matrix._ 
```C++
inline KOKKOS_FUNCTION DTensor < VectorIndexSet< Rho, Theta , Phi >, VectorIndexSet< R_cov , Z_cov, Zeta_cov > > ToroidalToCylindrical::inv_jacobian_matrix (
    CoordArg const & coord
) const
```



For different computations, we need either the complete Jacobian matrix or some of the coefficients of the matrix. The coefficients can be obtained independently with the function inv\_jacobian\_component.




**Parameters:**


* `coord` The coordinate where we evaluate the Jacobian matrix. 



**Returns:**

The Jacobian matrix. 





        

<hr>



### function jacobian 

_Compute the Jacobian: the determinant of the Jacobian matrix of the mapping._ 
```C++
inline KOKKOS_FUNCTION double ToroidalToCylindrical::jacobian (
    CoordArg const & coord
) const
```





**Parameters:**


* `coord` The coordinate where we evaluate the Jacobian.



**Returns:**

The value of the determinant of the Jacobian matrix. 





        

<hr>



### function jacobian\_component 

_Compute the (i,j) coefficient of the Jacobian matrix._ 
```C++
template<class IndexTag1, class IndexTag2>
inline KOKKOS_INLINE_FUNCTION double ToroidalToCylindrical::jacobian_component (
    CoordArg const & coord
) const
```





**Parameters:**


* `coord` The coordinate where we evaluate the Jacobian matrix.



**Template parameters:**


* `The` tag representing the index i. 
* `The` tag representing the index j.



**Returns:**

The value of the (i,j) coefficient of the Jacobian matrix. 





        

<hr>



### function jacobian\_matrix 

_Compute the Jacobian matrix._ 
```C++
inline KOKKOS_FUNCTION DTensor < VectorIndexSet< R , Z, Zeta >, VectorIndexSet< Rho_cov, Theta_cov , Phi_cov > > ToroidalToCylindrical::jacobian_matrix (
    CoordArg const & coord
) const
```



For different computations, we need either the complete Jacobian matrix or some of the coefficients of the matrix. The coefficients can be obtained independently with the function jacobian\_component.




**Parameters:**


* `coord` The coordinate where we evaluate the Jacobian matrix. 



**Returns:**

The Jacobian matrix. 





        

<hr>



### function operator() 

_Convert the_  _coordinate to the equivalent_ _coordinate._
```C++
inline KOKKOS_FUNCTION CoordResult ToroidalToCylindrical::operator() (
    CoordArg const & coord
) const
```





**Parameters:**


* `coord` The coordinate to be converted.



**Returns:**

The equivalent coordinate. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/mapping/toroidal_to_cylindrical.hpp`

