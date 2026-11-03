

# Struct connectivity\_details::FindInterface&lt; Edge, ddc::TypeSeq&lt; Interface1, RemainingInterfaceTypes... &gt; &gt;

**template &lt;class [**Edge**](structEdge.md), class Interface1, class... RemainingInterfaceTypes&gt;**



[**ClassList**](annotated.md) **>** [**connectivity\_details**](namespaceconnectivity__details.md) **>** [**FindInterface&lt; Edge, ddc::TypeSeq&lt; Interface1, RemainingInterfaceTypes... &gt; &gt;**](structconnectivity__details_1_1FindInterface_3_01Edge_00_01ddc_1_1TypeSeq_3_01Interface1_00_01Re9a613ed41c91d250bbda6f232579e462.md)



_Specialisation of_ [_**FindInterface**_](structconnectivity__details_1_1FindInterface.md) _to iterate recursively over the interface type sequence._

* `#include <connectivity_details.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef typename [**FindInterface**](structconnectivity__details_1_1FindInterface.md)&lt; [**Edge**](structEdge.md), ddc::TypeSeq&lt; RemainingInterfaceTypes... &gt; &gt;::type | [**type**](#typedef-type)  <br>_The type found by the class._  |
















































## Public Types Documentation




### typedef type 

_The type found by the class._ 
```C++
using connectivity_details::FindInterface< Edge, ddc::TypeSeq< Interface1, RemainingInterfaceTypes... > >::type =  typename FindInterface<Edge, ddc::TypeSeq<RemainingInterfaceTypes...>>::type;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/connectivity/connectivity_details.hpp`

