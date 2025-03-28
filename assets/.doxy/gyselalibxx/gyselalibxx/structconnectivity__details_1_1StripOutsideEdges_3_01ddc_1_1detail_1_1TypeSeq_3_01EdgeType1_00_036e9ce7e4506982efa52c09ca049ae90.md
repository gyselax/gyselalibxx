

# Struct connectivity\_details::StripOutsideEdges&lt; ddc::detail::TypeSeq&lt; EdgeType1, RemainingEdgeTypes... &gt; &gt;

**template &lt;class EdgeType1, class... RemainingEdgeTypes&gt;**



[**ClassList**](annotated.md) **>** [**connectivity\_details**](namespaceconnectivity__details.md) **>** [**StripOutsideEdges&lt; ddc::detail::TypeSeq&lt; EdgeType1, RemainingEdgeTypes... &gt; &gt;**](structconnectivity__details_1_1StripOutsideEdges_3_01ddc_1_1detail_1_1TypeSeq_3_01EdgeType1_00_036e9ce7e4506982efa52c09ca049ae90.md)



_Specialisation of_ [_**StripOutsideEdges**_](structconnectivity__details_1_1StripOutsideEdges.md) _to iterate recursively over the edge type sequence._

* `#include <connectivity_details.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef ddc::type\_seq\_merge\_t&lt; typename [**StripOutsideEdges**](structconnectivity__details_1_1StripOutsideEdges.md)&lt; ddc::detail::TypeSeq&lt; EdgeType1 &gt; &gt;::type, typename [**StripOutsideEdges**](structconnectivity__details_1_1StripOutsideEdges.md)&lt; ddc::detail::TypeSeq&lt; RemainingEdgeTypes... &gt; &gt;::type &gt; | [**type**](#typedef-type)  <br>_The type found by the class._  |
















































## Public Types Documentation




### typedef type 

_The type found by the class._ 
```C++
using connectivity_details::StripOutsideEdges< ddc::detail::TypeSeq< EdgeType1, RemainingEdgeTypes... > >::type =  ddc::type_seq_merge_t< typename StripOutsideEdges<ddc::detail::TypeSeq<EdgeType1> >::type, typename StripOutsideEdges<ddc::detail::TypeSeq<RemainingEdgeTypes...> >::type>;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/connectivity/connectivity_details.hpp`

