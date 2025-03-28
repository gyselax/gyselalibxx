

# Struct connectivity\_details::CollectInterfacesAlongDim&lt; StartEdge, InterfaceTypeSeq, insert\_pos, FoundInterfaces, MatchingEdge, false &gt;

**template &lt;class StartEdge, class InterfaceTypeSeq, InsertPosition insert\_pos, class FoundInterfaces, class MatchingEdge&gt;**



[**ClassList**](annotated.md) **>** [**connectivity\_details**](namespaceconnectivity__details.md) **>** [**CollectInterfacesAlongDim&lt; StartEdge, InterfaceTypeSeq, insert\_pos, FoundInterfaces, MatchingEdge, false &gt;**](structconnectivity__details_1_1CollectInterfacesAlongDim_3_01StartEdge_00_01InterfaceTypeSeq_00_b2108f65f3430e895714f416a2f43701.md)



_Specialisation of_ [_**CollectGridsAlongDim**_](structconnectivity__details_1_1CollectGridsAlongDim.md) _to iterate recursively over the grids on the dimension._

* `#include <connectivity_details.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef typename [**AddToTypeSeq**](structconnectivity__details_1_1AddToTypeSeq.md)&lt; enforce\_first\_interface\_edge\_t&lt; typename [**FindInterface**](structconnectivity__details_1_1FindInterface.md)&lt; StartEdge, InterfaceTypeSeq &gt;[**::type**](structconnectivity__details_1_1CollectInterfacesAlongDim_3_01StartEdge_00_01InterfaceTypeSeq_00_b2108f65f3430e895714f416a2f43701.md#typedef-type), std::conditional\_t&lt; StartEdge::extremity==FRONT, MatchingEdge, StartEdge &gt; &gt;, FoundInterfaces, insert\_pos &gt;[**::type**](structconnectivity__details_1_1CollectInterfacesAlongDim_3_01StartEdge_00_01InterfaceTypeSeq_00_b2108f65f3430e895714f416a2f43701.md#typedef-type) | [**NewInterfaceList**](#typedef-newinterfacelist)  <br>_The new list of interfaces that have been found including the interface from the current patch._  |
| typedef typename [**CollectInterfacesAlongDim**](structconnectivity__details_1_1CollectInterfacesAlongDim.md)&lt; typename [**SwapExtremity**](structconnectivity__details_1_1SwapExtremity.md)&lt; MatchingEdge &gt;::type, InterfaceTypeSeq, insert\_pos, [**NewInterfaceList**](structconnectivity__details_1_1CollectInterfacesAlongDim_3_01StartEdge_00_01InterfaceTypeSeq_00_b2108f65f3430e895714f416a2f43701.md#typedef-newinterfacelist) &gt;::type | [**type**](#typedef-type)  <br>_The type found by the class._  |
















































## Public Types Documentation




### typedef NewInterfaceList 

_The new list of interfaces that have been found including the interface from the current patch._ 
```C++
using connectivity_details::CollectInterfacesAlongDim< StartEdge, InterfaceTypeSeq, insert_pos, FoundInterfaces, MatchingEdge, false >::NewInterfaceList =  typename AddToTypeSeq< enforce_first_interface_edge_t< typename FindInterface<StartEdge, InterfaceTypeSeq>::type, std::conditional_t<StartEdge::extremity == FRONT, MatchingEdge, StartEdge> >, FoundInterfaces, insert_pos>::type;
```




<hr>



### typedef type 

_The type found by the class._ 
```C++
using connectivity_details::CollectInterfacesAlongDim< StartEdge, InterfaceTypeSeq, insert_pos, FoundInterfaces, MatchingEdge, false >::type =  typename CollectInterfacesAlongDim< typename SwapExtremity<MatchingEdge>::type, InterfaceTypeSeq, insert_pos, NewInterfaceList>::type;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/connectivity/connectivity_details.hpp`

