

# Struct connectivity\_details::StripOutsideEdges&lt; ddc::detail::TypeSeq&lt; EdgeType &gt; &gt;

**template &lt;class EdgeType&gt;**



[**ClassList**](annotated.md) **>** [**connectivity\_details**](namespaceconnectivity__details.md) **>** [**StripOutsideEdges&lt; ddc::detail::TypeSeq&lt; EdgeType &gt; &gt;**](structconnectivity__details_1_1StripOutsideEdges_3_01ddc_1_1detail_1_1TypeSeq_3_01EdgeType_01_4_01_4.md)



_Specialisation of_ [_**StripOutsideEdges**_](structconnectivity__details_1_1StripOutsideEdges.md) _for the case with one edge in the list._

* `#include <connectivity_details.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef std::conditional\_t&lt; std::is\_same\_v&lt; EdgeType, [**OutsideEdge**](structOutsideEdge.md) &gt;, ddc::detail::TypeSeq&lt;&gt;, ddc::detail::TypeSeq&lt; EdgeType &gt; &gt; | [**type**](#typedef-type)  <br>_The type found by the class._  |
















































## Public Types Documentation




### typedef type 

_The type found by the class._ 
```C++
using connectivity_details::StripOutsideEdges< ddc::detail::TypeSeq< EdgeType > >::type =  std::conditional_t< std::is_same_v<EdgeType, OutsideEdge>, ddc::detail::TypeSeq<>, ddc::detail::TypeSeq<EdgeType> >;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/connectivity/connectivity_details.hpp`

