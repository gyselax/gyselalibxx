

# Struct connectivity\_details::FindRelevantIdxRangeType&lt; QueryGrid1D, std::tuple&lt; IdxRangeHead, IdxRangeTypes... &gt; &gt;

**template &lt;class QueryGrid1D, class IdxRangeHead, class... IdxRangeTypes&gt;**



[**ClassList**](annotated.md) **>** [**connectivity\_details**](namespaceconnectivity__details.md) **>** [**FindRelevantIdxRangeType&lt; QueryGrid1D, std::tuple&lt; IdxRangeHead, IdxRangeTypes... &gt; &gt;**](structconnectivity__details_1_1FindRelevantIdxRangeType_3_01QueryGrid1D_00_01std_1_1tuple_3_01Id3b131c802b30082f4412eb4689d6d53b.md)



_Specialisation of_ [_**FindRelevantIdxRangeType**_](structconnectivity__details_1_1FindRelevantIdxRangeType.md) _to iterate recursively over the possible index range types._

* `#include <connectivity_details.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef ddc::type\_seq\_merge\_t&lt; typename [**SelectRelevantIdxRangeType**](structconnectivity__details_1_1SelectRelevantIdxRangeType.md)&lt; QueryGrid1D, IdxRangeHead &gt;::type, typename [**FindRelevantIdxRangeType**](structconnectivity__details_1_1FindRelevantIdxRangeType.md)&lt; QueryGrid1D, std::tuple&lt; IdxRangeTypes... &gt; &gt;::type &gt; | [**type**](#typedef-type)  <br>_The type found by the class._  |
















































## Public Types Documentation




### typedef type 

_The type found by the class._ 
```C++
using connectivity_details::FindRelevantIdxRangeType< QueryGrid1D, std::tuple< IdxRangeHead, IdxRangeTypes... > >::type =  ddc::type_seq_merge_t< typename SelectRelevantIdxRangeType<QueryGrid1D, IdxRangeHead>::type, typename FindRelevantIdxRangeType<QueryGrid1D, std::tuple<IdxRangeTypes...> >::type>;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/connectivity/connectivity_details.hpp`

