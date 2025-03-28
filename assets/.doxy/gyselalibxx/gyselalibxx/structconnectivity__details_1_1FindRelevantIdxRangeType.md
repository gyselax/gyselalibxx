

# Struct connectivity\_details::FindRelevantIdxRangeType

**template &lt;class QueryGrid1D, class IdxRangeTuple&gt;**



[**ClassList**](annotated.md) **>** [**connectivity\_details**](namespaceconnectivity__details.md) **>** [**FindRelevantIdxRangeType**](structconnectivity__details_1_1FindRelevantIdxRangeType.md)



_A class to find any index range types which contain an index range defined on the provided grid. E.g. Grid1, std::tuple&lt;IdxRange&lt;Grid1, Grid2&gt;, IdxRange&lt;Grid3,Grid4&gt;&gt; will find: ddc::detail::TypeSeq&lt;IdxRange&lt;Grid1, Grid2&gt;&gt;_ [More...](#detailed-description)


































































## Detailed Description




**Template parameters:**


* `QueryGrid1D` The grid being searched for. 
* `A` tuple of index ranges which may contain the relevant grid. 




    

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/connectivity/connectivity_details.hpp`

