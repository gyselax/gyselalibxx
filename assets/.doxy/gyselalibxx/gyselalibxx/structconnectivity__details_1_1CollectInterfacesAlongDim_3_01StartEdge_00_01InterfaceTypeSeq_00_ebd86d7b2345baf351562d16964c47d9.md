

# Struct connectivity\_details::CollectInterfacesAlongDim&lt; StartEdge, InterfaceTypeSeq, insert\_pos, FoundInterfaces, OutsideEdge, false &gt;

**template &lt;class StartEdge, class InterfaceTypeSeq, InsertPosition insert\_pos, class FoundInterfaces&gt;**



[**ClassList**](annotated.md) **>** [**connectivity\_details**](namespaceconnectivity__details.md) **>** [**CollectInterfacesAlongDim&lt; StartEdge, InterfaceTypeSeq, insert\_pos, FoundInterfaces, OutsideEdge, false &gt;**](structconnectivity__details_1_1CollectInterfacesAlongDim_3_01StartEdge_00_01InterfaceTypeSeq_00_ebd86d7b2345baf351562d16964c47d9.md)



_Specialisation of_ [_**CollectGridsAlongDim**_](structconnectivity__details_1_1CollectGridsAlongDim.md) _to stop when there are no more grids._

* `#include <connectivity_details.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef typename [**AddToTypeSeq**](structconnectivity__details_1_1AddToTypeSeq.md)&lt; enforce\_first\_interface\_edge\_t&lt; typename [**FindInterface**](structconnectivity__details_1_1FindInterface.md)&lt; StartEdge, InterfaceTypeSeq &gt;::type, std::conditional\_t&lt; StartEdge::extremity==FRONT, [**OutsideEdge**](structOutsideEdge.md), StartEdge &gt; &gt;, FoundInterfaces, insert\_pos &gt;::type | [**type**](#typedef-type)  <br>_The type found by the class._  |
















































## Public Types Documentation




### typedef type 

_The type found by the class._ 
```C++
using connectivity_details::CollectInterfacesAlongDim< StartEdge, InterfaceTypeSeq, insert_pos, FoundInterfaces, OutsideEdge, false >::type =  typename AddToTypeSeq< enforce_first_interface_edge_t< typename FindInterface<StartEdge, InterfaceTypeSeq>::type, std::conditional_t<StartEdge::extremity == FRONT, OutsideEdge, StartEdge> >, FoundInterfaces, insert_pos>::type;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/connectivity/connectivity_details.hpp`

