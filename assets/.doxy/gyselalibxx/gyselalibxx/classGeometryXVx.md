

# Class GeometryXVx



[**ClassList**](annotated.md) **>** [**GeometryXVx**](classGeometryXVx.md)



_A class providing aliases for useful subindex ranges of the geometry. It is used as template parameter for generic dimensionality-agnostic operators such as advections._ 

* `#include <geometry.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef IdxRangeSpXVx | [**IdxRangeFdistribu**](#typedef-idxrangefdistribu)  <br>_An alias for the whole distribution function discrete index range type._  |
| typedef IdxRangeX | [**IdxRangeSpatial**](#typedef-idxrangespatial)  <br>_An alias for the spatial discrete index range type._  |
| typedef IdxRangeVx | [**IdxRangeVelocity**](#typedef-idxrangevelocity)  <br>_An alias for the velocity discrete index range type._  |
| typedef std::conditional\_t&lt; std::is\_same\_v&lt; [**T**](structT.md), [**GridVx**](structGridVx.md) &gt;, [**GridX**](structGridX.md), void &gt; | [**spatial\_dim\_for**](#typedef-spatial_dim_for)  <br>_A templated type giving the spatial discretised dimension type associated to a velocity discretised dimension type._  |
| typedef std::conditional\_t&lt; std::is\_same\_v&lt; [**T**](structT.md), [**GridX**](structGridX.md) &gt;, [**GridVx**](structGridVx.md), void &gt; | [**velocity\_dim\_for**](#typedef-velocity_dim_for)  <br>_A templated type giving the velocity discretised dimension type associated to a spatial discretised dimension type._  |
















































## Public Types Documentation




### typedef IdxRangeFdistribu 

_An alias for the whole distribution function discrete index range type._ 
```C++
using GeometryXVx::IdxRangeFdistribu =  IdxRangeSpXVx;
```




<hr>



### typedef IdxRangeSpatial 

_An alias for the spatial discrete index range type._ 
```C++
using GeometryXVx::IdxRangeSpatial =  IdxRangeX;
```




<hr>



### typedef IdxRangeVelocity 

_An alias for the velocity discrete index range type._ 
```C++
using GeometryXVx::IdxRangeVelocity =  IdxRangeVx;
```




<hr>



### typedef spatial\_dim\_for 

_A templated type giving the spatial discretised dimension type associated to a velocity discretised dimension type._ 
```C++
using GeometryXVx::spatial_dim_for =  std::conditional_t<std::is_same_v<T, GridVx>, GridX, void>;
```




<hr>



### typedef velocity\_dim\_for 

_A templated type giving the velocity discretised dimension type associated to a spatial discretised dimension type._ 
```C++
using GeometryXVx::velocity_dim_for =  std::conditional_t<std::is_same_v<T, GridX>, GridVx, void>;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/geometryXVx/geometry/geometry.hpp`

