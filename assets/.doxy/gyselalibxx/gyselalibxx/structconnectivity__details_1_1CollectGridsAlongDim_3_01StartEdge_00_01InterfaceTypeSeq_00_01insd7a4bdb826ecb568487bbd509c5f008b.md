

# Struct connectivity\_details::CollectGridsAlongDim&lt; StartEdge, InterfaceTypeSeq, insert\_pos, FoundGrids, MatchingEdge, false &gt;

**template &lt;class StartEdge, class InterfaceTypeSeq, InsertPosition insert\_pos, class FoundGrids, class MatchingEdge&gt;**



[**ClassList**](annotated.md) **>** [**connectivity\_details**](namespaceconnectivity__details.md) **>** [**CollectGridsAlongDim&lt; StartEdge, InterfaceTypeSeq, insert\_pos, FoundGrids, MatchingEdge, false &gt;**](structconnectivity__details_1_1CollectGridsAlongDim_3_01StartEdge_00_01InterfaceTypeSeq_00_01insd7a4bdb826ecb568487bbd509c5f008b.md)



_Specialisation of_ [_**CollectGridsAlongDim**_](structconnectivity__details_1_1CollectGridsAlongDim.md) _to iterate recursively over the grids on the dimension._

* `#include <connectivity_details.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef typename [**AddToTypeSeq**](structconnectivity__details_1_1AddToTypeSeq.md)&lt; typename StartEdge::perpendicular\_grid, FoundGrids, insert\_pos &gt;[**::type**](structconnectivity__details_1_1CollectGridsAlongDim_3_01StartEdge_00_01InterfaceTypeSeq_00_01insd7a4bdb826ecb568487bbd509c5f008b.md#typedef-type) | [**NewGridList**](#typedef-newgridlist)  <br>_The new list of grids that have been found including the grid from the current patch._  |
| typedef typename [**CollectGridsAlongDim**](structconnectivity__details_1_1CollectGridsAlongDim.md)&lt; typename [**SwapExtremity**](structconnectivity__details_1_1SwapExtremity.md)&lt; MatchingEdge &gt;::type, InterfaceTypeSeq, insert\_pos, [**NewGridList**](structconnectivity__details_1_1CollectGridsAlongDim_3_01StartEdge_00_01InterfaceTypeSeq_00_01insd7a4bdb826ecb568487bbd509c5f008b.md#typedef-newgridlist) &gt;::type | [**type**](#typedef-type)  <br>_The type found by the class._  |
















































## Public Types Documentation




### typedef NewGridList 

_The new list of grids that have been found including the grid from the current patch._ 
```C++
using connectivity_details::CollectGridsAlongDim< StartEdge, InterfaceTypeSeq, insert_pos, FoundGrids, MatchingEdge, false >::NewGridList =  typename AddToTypeSeq<typename StartEdge::perpendicular_grid, FoundGrids, insert_pos>:: type;
```




<hr>



### typedef type 

_The type found by the class._ 
```C++
using connectivity_details::CollectGridsAlongDim< StartEdge, InterfaceTypeSeq, insert_pos, FoundGrids, MatchingEdge, false >::type =  typename CollectGridsAlongDim< typename SwapExtremity<MatchingEdge>::type, InterfaceTypeSeq, insert_pos, NewGridList>::type;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/connectivity/connectivity_details.hpp`

