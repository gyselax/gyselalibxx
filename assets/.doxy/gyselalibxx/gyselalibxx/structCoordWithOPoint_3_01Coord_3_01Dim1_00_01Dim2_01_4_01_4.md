

# Struct CoordWithOPoint&lt; Coord&lt; Dim1, Dim2 &gt; &gt;

**template &lt;class Dim1, class Dim2&gt;**



[**ClassList**](annotated.md) **>** [**CoordWithOPoint&lt; Coord&lt; Dim1, Dim2 &gt; &gt;**](structCoordWithOPoint_3_01Coord_3_01Dim1_00_01Dim2_01_4_01_4.md)



_A class to identify the radial and poloidal components of a 2D cuvilinear coordinate._ 

* `#include <coord_transformation_tools.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef std::conditional\_t&lt; Dim1::PERIODIC, Dim2, Dim1 &gt; | [**curvilinear\_tag\_r**](#typedef-curvilinear_tag_r)  <br>_Radial tag._  |
| typedef std::conditional\_t&lt; Dim1::PERIODIC, Dim1, Dim2 &gt; | [**curvilinear\_tag\_theta**](#typedef-curvilinear_tag_theta)  <br>_Poloidal tag._  |
















































## Public Types Documentation




### typedef curvilinear\_tag\_r 

_Radial tag._ 
```C++
using CoordWithOPoint< Coord< Dim1, Dim2 > >::curvilinear_tag_r =  std::conditional_t<Dim1::PERIODIC, Dim2, Dim1>;
```




<hr>



### typedef curvilinear\_tag\_theta 

_Poloidal tag._ 
```C++
using CoordWithOPoint< Coord< Dim1, Dim2 > >::curvilinear_tag_theta =  std::conditional_t<Dim1::PERIODIC, Dim1, Dim2>;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/coord_transformations/coord_transformation_tools.hpp`

