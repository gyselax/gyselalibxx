

# Class GMGPolarTools::MappingToDomainGeometry

**template &lt;class ToPhysicalMapping&gt;**



[**ClassList**](annotated.md) **>** [**GMGPolarTools**](namespaceGMGPolarTools.md) **>** [**MappingToDomainGeometry**](classGMGPolarTools_1_1MappingToDomainGeometry.md)



_Wraps a gyselalibxx coordinate mapping to satisfy the GMGPolar DomainGeometry concept._ [More...](#detailed-description)

* `#include <gmg_polar_poisson_like_solver.hpp>`





































## Public Functions

| Type | Name |
| ---: | :--- |
|  KOKKOS\_INLINE\_FUNCTION double | [**Fx**](#function-fx) (const double & r, const double & theta) const<br>_X(r, theta)_  |
|  KOKKOS\_INLINE\_FUNCTION double | [**Fy**](#function-fy) (const double & r, const double & theta) const<br>_Y(r, theta)_  |
|   | [**MappingToDomainGeometry**](#function-mappingtodomaingeometry) (ToPhysicalMapping to\_physical) <br>_Construct the wrapper class._  |
|  KOKKOS\_INLINE\_FUNCTION double | [**dFx\_dr**](#function-dfx_dr) (const double & r, const double & theta) const<br>_d/dr X(r, theta)_  |
|  KOKKOS\_INLINE\_FUNCTION double | [**dFx\_dtheta**](#function-dfx_dtheta) (const double & r, const double & theta) const<br>_d/(d theta) X(r, theta)_  |
|  KOKKOS\_INLINE\_FUNCTION double | [**dFy\_dr**](#function-dfy_dr) (const double & r, const double & theta) const<br>_d/dr Y(r, theta)_  |
|  KOKKOS\_INLINE\_FUNCTION double | [**dFy\_dtheta**](#function-dfy_dtheta) (const double & r, const double & theta) const<br>_d/(d theta) Y(r, theta)_  |




























## Detailed Description




**Template parameters:**


* `ToPhysicalMapping` A mapping from (r, theta) curvilinear coordinates to (x, y) Cartesian. 




    
## Public Functions Documentation




### function Fx 

_X(r, theta)_ 
```C++
inline KOKKOS_INLINE_FUNCTION double GMGPolarTools::MappingToDomainGeometry::Fx (
    const double & r,
    const double & theta
) const
```




<hr>



### function Fy 

_Y(r, theta)_ 
```C++
inline KOKKOS_INLINE_FUNCTION double GMGPolarTools::MappingToDomainGeometry::Fy (
    const double & r,
    const double & theta
) const
```




<hr>



### function MappingToDomainGeometry 

_Construct the wrapper class._ 
```C++
inline explicit GMGPolarTools::MappingToDomainGeometry::MappingToDomainGeometry (
    ToPhysicalMapping to_physical
) 
```




<hr>



### function dFx\_dr 

_d/dr X(r, theta)_ 
```C++
inline KOKKOS_INLINE_FUNCTION double GMGPolarTools::MappingToDomainGeometry::dFx_dr (
    const double & r,
    const double & theta
) const
```




<hr>



### function dFx\_dtheta 

_d/(d theta) X(r, theta)_ 
```C++
inline KOKKOS_INLINE_FUNCTION double GMGPolarTools::MappingToDomainGeometry::dFx_dtheta (
    const double & r,
    const double & theta
) const
```




<hr>



### function dFy\_dr 

_d/dr Y(r, theta)_ 
```C++
inline KOKKOS_INLINE_FUNCTION double GMGPolarTools::MappingToDomainGeometry::dFy_dr (
    const double & r,
    const double & theta
) const
```




<hr>



### function dFy\_dtheta 

_d/(d theta) Y(r, theta)_ 
```C++
inline KOKKOS_INLINE_FUNCTION double GMGPolarTools::MappingToDomainGeometry::dFy_dtheta (
    const double & r,
    const double & theta
) const
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/pde_solvers/gmg_polar_poisson_like_solver.hpp`

