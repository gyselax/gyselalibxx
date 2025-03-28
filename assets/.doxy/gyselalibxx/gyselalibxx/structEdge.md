

# Struct Edge

**template &lt;class [**Patch**](structPatch.md), class Grid1D, Extremity extremity\_val&gt;**



[**ClassList**](annotated.md) **>** [**Edge**](structEdge.md)



_Define an edge of a given patch._ [More...](#detailed-description)

* `#include <edge.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef [**Patch**](structPatch.md) | [**associated\_patch**](#typedef-associated_patch)  <br>[_**Patch**_](structPatch.md) _where the edge is defined._ |
| typedef std::conditional\_t&lt; std::is\_same\_v&lt; Grid1D, typename Patch::Grid1 &gt;, typename Patch::Grid2, typename Patch::Grid1 &gt; | [**parallel\_grid**](#typedef-parallel_grid)  <br>_Grid parallel to the edge._  |
| typedef Grid1D | [**perpendicular\_grid**](#typedef-perpendicular_grid)  <br>_Grid on the perpendicular dimension of the edge._  |






## Public Static Attributes

| Type | Name |
| ---: | :--- |
|  constexpr Extremity | [**extremity**](#variable-extremity)   = `extremity\_val`<br>_Design if the edge is on the BACK or the FRONT of the other dimension._  |










































## Detailed Description


An edge is defined by a patch, a dimension and an extremity. For example, in the patch defined on logical index range ,



* the edge [**GridY**](structGridY.md), BACK refers to the set ,
* and the edge [**GridY**](structGridY.md), FRONT refers to the set .






**Template parameters:**


* [**Patch**](structPatch.md) [**Patch**](structPatch.md) where the edge is defined. 
* `Grid1D` Grid on the complementary dimension of the edge. 
 
* `extremity_val` The BACK or FRONT value. 




    
## Public Types Documentation




### typedef associated\_patch 

[_**Patch**_](structPatch.md) _where the edge is defined._
```C++
using Edge< Patch, Grid1D, extremity_val >::associated_patch =  Patch;
```




<hr>



### typedef parallel\_grid 

_Grid parallel to the edge._ 
```C++
using Edge< Patch, Grid1D, extremity_val >::parallel_grid =  std::conditional_t< std::is_same_v<Grid1D, typename Patch::Grid1>, typename Patch::Grid2, typename Patch::Grid1>;
```




<hr>



### typedef perpendicular\_grid 

_Grid on the perpendicular dimension of the edge._ 
```C++
using Edge< Patch, Grid1D, extremity_val >::perpendicular_grid =  Grid1D;
```




<hr>
## Public Static Attributes Documentation




### variable extremity 

_Design if the edge is on the BACK or the FRONT of the other dimension._ 
```C++
constexpr Extremity Edge< Patch, Grid1D, extremity_val >::extremity;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/connectivity/edge.hpp`

