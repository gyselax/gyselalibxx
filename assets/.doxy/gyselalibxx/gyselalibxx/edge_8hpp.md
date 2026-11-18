

# File edge.hpp



[**FileList**](files.md) **>** [**connectivity**](dir_28b51abc9241105ab41b66c468e7d019.md) **>** [**edge.hpp**](edge_8hpp.md)

[Go to the source code of this file](edge_8hpp_source.md)



* `#include <ddc/ddc.hpp>`
* `#include "ddc_aliases.hpp"`
* `#include "geometry_descriptors.hpp"`
* `#include "patch.hpp"`













## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**concepts**](namespaceconcepts.md) <br> |


## Classes

| Type | Name |
| ---: | :--- |
| struct | [**Edge**](structEdge.md) &lt;Patch, class Grid1D, extremity\_val&gt;<br>_Define an edge of a given patch._  |
| struct | [**OutsideEdge**](structOutsideEdge.md) <br>_Define an edge for the outside index range._ [_**OutsideEdge**_](structOutsideEdge.md) _is a pseudo-edge outside the index range used to define interfaces between patches and the outside index range._ |
| struct | [**is\_edge**](structis__edge.md) &lt;class EdgeType&gt;<br>_A type trait to check if a type is an edge (an_ [_**Edge**_](structEdge.md) _or an_[_**OutsideEdge**_](structOutsideEdge.md) _)._ |
| struct | [**is\_edge&lt; Edge&lt; Patch, Grid1D, extremity\_val &gt; &gt;**](structis__edge_3_01Edge_3_01Patch_00_01Grid1D_00_01extremity__val_01_4_01_4.md) &lt;Patch, class Grid1D, extremity\_val&gt;<br>_Specialisation of_ [_**is\_edge**_](structis__edge.md) _for an_[_**Edge**_](structEdge.md) _of a patch._ |
| struct | [**is\_edge&lt; OutsideEdge &gt;**](structis__edge_3_01OutsideEdge_01_4.md) &lt;&gt;<br>_Specialisation of_ [_**is\_edge**_](structis__edge.md) _for the_[_**OutsideEdge**_](structOutsideEdge.md) _._ |






## Public Attributes

| Type | Name |
| ---: | :--- |
|  bool constexpr | [**is\_edge\_v**](#variable-is_edge_v)   = `[**is\_edge**](structis__edge.md)&lt;EdgeType&gt;::value`<br>_A boolean indicating if a type is an edge (an_ [_**Edge**_](structEdge.md) _or an_[_**OutsideEdge**_](structOutsideEdge.md) _)._ |
















## Public Functions

| Type | Name |
| ---: | :--- |
|  int | [**get\_edge\_patch\_id**](#function-get_edge_patch_id) () <br>_Get the id of the patch on which an edge is defined._  |




























## Public Attributes Documentation




### variable is\_edge\_v 

_A boolean indicating if a type is an edge (an_ [_**Edge**_](structEdge.md) _or an_[_**OutsideEdge**_](structOutsideEdge.md) _)._
```C++
bool constexpr is_edge_v;
```





**Template parameters:**


* `EdgeType` The type being checked. 




        

<hr>
## Public Functions Documentation




### function get\_edge\_patch\_id 

_Get the id of the patch on which an edge is defined._ 
```C++
template<class EdgeType>
int get_edge_patch_id () 
```



This can be used to obtain a compile-time ordering of the edges.




**Template parameters:**


* `EdgeType` The edge (an [**Edge**](structEdge.md) or an [**OutsideEdge**](structOutsideEdge.md)). 



**Returns:**

The id of the patch associated with the edge, or -1 for an [**OutsideEdge**](structOutsideEdge.md). 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/connectivity/edge.hpp`

