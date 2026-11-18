

# Struct connectivity\_details::FindInterface&lt; Edge, ddc::TypeSeq&lt;::Interface&lt; OEdge, Edge, Orientations &gt;, RemainingInterfaceTypes... &gt; &gt;

**template &lt;concepts::Edge Edge, concepts::Edge OEdge, bool Orientations, class... RemainingInterfaceTypes&gt;**



[**ClassList**](annotated.md) **>** [**connectivity\_details**](namespaceconnectivity__details.md) **>** [**FindInterface&lt; Edge, ddc::TypeSeq&lt;::Interface&lt; OEdge, Edge, Orientations &gt;, RemainingInterfaceTypes... &gt; &gt;**](structconnectivity__details_1_1FindInterface_3_01Edge_00_01ddc_1_1TypeSeq_3_1_1Interface_3_01OEdcf579b47e460dc35ee239300e37577bd.md)



_Specialisation of_ [_**FindInterface**_](structconnectivity__details_1_1FindInterface.md) _for the case where Edge1 from the second interface matches_[_**Edge**_](structEdge.md) _._

* `#include <connectivity_details.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef [**Interface**](structInterface.md)&lt; OEdge, [**Edge**](structEdge.md), Orientations &gt; | [**type**](#typedef-type)  <br>_The type found by the class._  |
















































## Public Types Documentation




### typedef type 

_The type found by the class._ 
```C++
using connectivity_details::FindInterface< Edge, ddc::TypeSeq<::Interface< OEdge, Edge, Orientations >, RemainingInterfaceTypes... > >::type =  Interface<OEdge, Edge, Orientations>;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/connectivity/connectivity_details.hpp`

