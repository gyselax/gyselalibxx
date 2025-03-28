

# Struct connectivity\_details::ExtractPatches&lt; ddc::detail::TypeSeq&lt; EdgeType1, EdgeTypes... &gt; &gt;

**template &lt;class EdgeType1, class... EdgeTypes&gt;**



[**ClassList**](annotated.md) **>** [**connectivity\_details**](namespaceconnectivity__details.md) **>** [**ExtractPatches&lt; ddc::detail::TypeSeq&lt; EdgeType1, EdgeTypes... &gt; &gt;**](structconnectivity__details_1_1ExtractPatches_3_01ddc_1_1detail_1_1TypeSeq_3_01EdgeType1_00_01EdgeTypes_8_8_8_01_4_01_4.md)



_Specialisation of_ [_**ExtractPatches**_](structconnectivity__details_1_1ExtractPatches.md) _to iterate recursively over the edge type sequence._

* `#include <connectivity_details.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef ddc::type\_seq\_merge\_t&lt; ddc::detail::TypeSeq&lt; typename EdgeType1::associated\_patch &gt;, typename [**ExtractPatches**](structconnectivity__details_1_1ExtractPatches.md)&lt; ddc::detail::TypeSeq&lt; EdgeTypes... &gt; &gt;::type &gt; | [**type**](#typedef-type)  <br>_The type found by the class._  |
















































## Public Types Documentation




### typedef type 

_The type found by the class._ 
```C++
using connectivity_details::ExtractPatches< ddc::detail::TypeSeq< EdgeType1, EdgeTypes... > >::type =  ddc::type_seq_merge_t< ddc::detail::TypeSeq<typename EdgeType1::associated_patch>, typename ExtractPatches<ddc::detail::TypeSeq<EdgeTypes...> >::type>;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/connectivity/connectivity_details.hpp`

