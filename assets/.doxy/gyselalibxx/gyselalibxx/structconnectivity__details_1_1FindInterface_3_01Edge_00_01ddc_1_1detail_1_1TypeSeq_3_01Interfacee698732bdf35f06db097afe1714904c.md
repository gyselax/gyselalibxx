

# Struct connectivity\_details::FindInterface&lt; Edge, ddc::detail::TypeSeq&lt; Interface&lt; OEdge, Edge, Orientations &gt;, RemainingInterfaceTypes... &gt; &gt;

**template &lt;class [**Edge**](structEdge.md), class OEdge, bool Orientations, class... RemainingInterfaceTypes&gt;**



[**ClassList**](annotated.md) **>** [**connectivity\_details**](namespaceconnectivity__details.md) **>** [**FindInterface&lt; Edge, ddc::detail::TypeSeq&lt; Interface&lt; OEdge, Edge, Orientations &gt;, RemainingInterfaceTypes... &gt; &gt;**](structconnectivity__details_1_1FindInterface_3_01Edge_00_01ddc_1_1detail_1_1TypeSeq_3_01Interfacee698732bdf35f06db097afe1714904c.md)



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
using connectivity_details::FindInterface< Edge, ddc::detail::TypeSeq< Interface< OEdge, Edge, Orientations >, RemainingInterfaceTypes... > >::type =  Interface<OEdge, Edge, Orientations>;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/connectivity/connectivity_details.hpp`

