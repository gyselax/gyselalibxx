

# Struct connectivity\_details::SelectRelevantIdxRangeType&lt; QueryGrid1D, IdxRange&lt; IdxRangeGrids... &gt; &gt;

**template &lt;class QueryGrid1D, class... IdxRangeGrids&gt;**



[**ClassList**](annotated.md) **>** [**connectivity\_details**](namespaceconnectivity__details.md) **>** [**SelectRelevantIdxRangeType&lt; QueryGrid1D, IdxRange&lt; IdxRangeGrids... &gt; &gt;**](structconnectivity__details_1_1SelectRelevantIdxRangeType_3_01QueryGrid1D_00_01IdxRange_3_01IdxRangeGrids_8_8_8_01_4_01_4.md)



_Specialisation of_ [_**SelectRelevantIdxRangeType**_](structconnectivity__details_1_1SelectRelevantIdxRangeType.md) _to get access to the grids in the index range._

* `#include <connectivity_details.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef std::conditional\_t&lt; ddc::in\_tags\_v&lt; QueryGrid1D, ddc::detail::TypeSeq&lt; IdxRangeGrids... &gt; &gt;, ddc::detail::TypeSeq&lt; IdxRange&lt; IdxRangeGrids... &gt; &gt;, ddc::detail::TypeSeq&lt;&gt; &gt; | [**type**](#typedef-type)  <br>_The type found by the class._  |
















































## Public Types Documentation




### typedef type 

_The type found by the class._ 
```C++
using connectivity_details::SelectRelevantIdxRangeType< QueryGrid1D, IdxRange< IdxRangeGrids... > >::type =  std::conditional_t< ddc::in_tags_v<QueryGrid1D, ddc::detail::TypeSeq<IdxRangeGrids...> >, ddc::detail::TypeSeq<IdxRange<IdxRangeGrids...> >, ddc::detail::TypeSeq<> >;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/connectivity/connectivity_details.hpp`

