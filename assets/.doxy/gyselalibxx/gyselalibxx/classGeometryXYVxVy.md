

# Class GeometryXYVxVy



[**ClassList**](annotated.md) **>** [**GeometryXYVxVy**](classGeometryXYVxVy.md)



_A class providing aliases for useful subindex ranges of the geometry when the data is saved with the spatial dimensions distributed across MPI ranks. It is used as template parameter for generic dimensionality-agnostic operators such as advections._ 

* `#include <geometry.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef IdxRangeSpXYVxVy | [**IdxRangeFdistribu**](#typedef-idxrangefdistribu)  <br>_An alias for the whole distribution function discrete index range type._  |
| typedef IdxRangeXY | [**IdxRangeSpatial**](#typedef-idxrangespatial)  <br>_A templated type giving the spatial discretised dimension type associated to a velocity discretised dimension type._  |
| typedef IdxRangeVxVy | [**IdxRangeVelocity**](#typedef-idxrangevelocity)  <br>_An alias for the velocity discrete index range type._  |
| typedef std::conditional\_t&lt; std::is\_same\_v&lt; [**T**](structT.md), [**GridX**](structGridX.md) &gt;, [**GridVx**](structGridVx.md), std::conditional\_t&lt; std::is\_same\_v&lt; [**T**](structT.md), [**GridY**](structGridY.md) &gt;, [**GridVy**](structGridVy.md), void &gt; &gt; | [**velocity\_dim\_for**](#typedef-velocity_dim_for)  <br>_A templated type giving the velocity discretised dimension type associated to a spatial discretised dimension type._  |
















































## Public Types Documentation




### typedef IdxRangeFdistribu 

_An alias for the whole distribution function discrete index range type._ 
```C++
using GeometryXYVxVy::IdxRangeFdistribu =  IdxRangeSpXYVxVy;
```




<hr>



### typedef IdxRangeSpatial 

_A templated type giving the spatial discretised dimension type associated to a velocity discretised dimension type._ 
```C++
using GeometryXYVxVy::IdxRangeSpatial =  IdxRangeXY;
```



An alias for the spatial discrete index range type. 


        

<hr>



### typedef IdxRangeVelocity 

_An alias for the velocity discrete index range type._ 
```C++
using GeometryXYVxVy::IdxRangeVelocity =  IdxRangeVxVy;
```




<hr>



### typedef velocity\_dim\_for 

_A templated type giving the velocity discretised dimension type associated to a spatial discretised dimension type._ 
```C++
using GeometryXYVxVy::velocity_dim_for =  std::conditional_t< std::is_same_v<T, GridX>, GridVx, std::conditional_t<std::is_same_v<T, GridY>, GridVy, void> >;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/geometryXYVxVy/geometry/geometry.hpp`

