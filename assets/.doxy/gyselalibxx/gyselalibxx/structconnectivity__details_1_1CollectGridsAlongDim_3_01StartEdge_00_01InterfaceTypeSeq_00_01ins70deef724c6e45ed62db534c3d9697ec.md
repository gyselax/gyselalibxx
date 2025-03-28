

# Struct connectivity\_details::CollectGridsAlongDim&lt; StartEdge, InterfaceTypeSeq, insert\_pos, FoundGrids, OutsideEdge, false &gt;

**template &lt;class StartEdge, class InterfaceTypeSeq, InsertPosition insert\_pos, class FoundGrids&gt;**



[**ClassList**](annotated.md) **>** [**connectivity\_details**](namespaceconnectivity__details.md) **>** [**CollectGridsAlongDim&lt; StartEdge, InterfaceTypeSeq, insert\_pos, FoundGrids, OutsideEdge, false &gt;**](structconnectivity__details_1_1CollectGridsAlongDim_3_01StartEdge_00_01InterfaceTypeSeq_00_01ins70deef724c6e45ed62db534c3d9697ec.md)



_Specialisation of_ [_**CollectGridsAlongDim**_](structconnectivity__details_1_1CollectGridsAlongDim.md) _to stop when there are no more grids._

* `#include <connectivity_details.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef typename [**AddToTypeSeq**](structconnectivity__details_1_1AddToTypeSeq.md)&lt; typename StartEdge::perpendicular\_grid, FoundGrids, insert\_pos &gt;::type | [**type**](#typedef-type)  <br>_The type found by the class._  |
















































## Public Types Documentation




### typedef type 

_The type found by the class._ 
```C++
using connectivity_details::CollectGridsAlongDim< StartEdge, InterfaceTypeSeq, insert_pos, FoundGrids, OutsideEdge, false >::type =  typename AddToTypeSeq<typename StartEdge::perpendicular_grid, FoundGrids, insert_pos>:: type;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/connectivity/connectivity_details.hpp`

