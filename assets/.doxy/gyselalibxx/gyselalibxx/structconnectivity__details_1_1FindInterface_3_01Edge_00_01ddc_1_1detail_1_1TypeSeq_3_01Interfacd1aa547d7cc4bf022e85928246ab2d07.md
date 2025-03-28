

# Struct connectivity\_details::FindInterface&lt; Edge, ddc::detail::TypeSeq&lt; Interface&lt; Edge, OEdge, Orientations &gt;, RemainingInterfaceTypes... &gt; &gt;

**template &lt;class [**Edge**](structEdge.md), class OEdge, bool Orientations, class... RemainingInterfaceTypes&gt;**



[**ClassList**](annotated.md) **>** [**connectivity\_details**](namespaceconnectivity__details.md) **>** [**FindInterface&lt; Edge, ddc::detail::TypeSeq&lt; Interface&lt; Edge, OEdge, Orientations &gt;, RemainingInterfaceTypes... &gt; &gt;**](structconnectivity__details_1_1FindInterface_3_01Edge_00_01ddc_1_1detail_1_1TypeSeq_3_01Interfacd1aa547d7cc4bf022e85928246ab2d07.md)



_Specialisation of_ [_**FindInterface**_](structconnectivity__details_1_1FindInterface.md) _for the case where Edge1 from the first interface matches_[_**Edge**_](structEdge.md) _._

* `#include <connectivity_details.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef [**Interface**](structInterface.md)&lt; [**Edge**](structEdge.md), OEdge, Orientations &gt; | [**type**](#typedef-type)  <br>_The type found by the class._  |
















































## Public Types Documentation




### typedef type 

_The type found by the class._ 
```C++
using connectivity_details::FindInterface< Edge, ddc::detail::TypeSeq< Interface< Edge, OEdge, Orientations >, RemainingInterfaceTypes... > >::type =  Interface<Edge, OEdge, Orientations>;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/connectivity/connectivity_details.hpp`

