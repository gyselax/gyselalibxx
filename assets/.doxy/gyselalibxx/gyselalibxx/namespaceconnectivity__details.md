

# Namespace connectivity\_details



[**Namespace List**](namespaces.md) **>** [**connectivity\_details**](namespaceconnectivity__details.md)




















## Classes

| Type | Name |
| ---: | :--- |
| struct | [**AddToTypeSeq**](structconnectivity__details_1_1AddToTypeSeq.md) &lt;class ToInsert, class TypeSeq, insert\_pos&gt;<br>_A class which helps insert an element into a type sequence._  |
| struct | [**AddToTypeSeq&lt; ToInsert, TypeSeq, BackInsert &gt;**](structconnectivity__details_1_1AddToTypeSeq_3_01ToInsert_00_01TypeSeq_00_01BackInsert_01_4.md) &lt;class ToInsert, class TypeSeq&gt;<br>_Specialisation of_ [_**AddToTypeSeq**_](structconnectivity__details_1_1AddToTypeSeq.md) _to add an element at the back of the type sequence._ |
| struct | [**AddToTypeSeq&lt; ToInsert, TypeSeq, FrontInsert &gt;**](structconnectivity__details_1_1AddToTypeSeq_3_01ToInsert_00_01TypeSeq_00_01FrontInsert_01_4.md) &lt;class ToInsert, class TypeSeq&gt;<br>_Specialisation of_ [_**AddToTypeSeq**_](structconnectivity__details_1_1AddToTypeSeq.md) _to add an element at the front of the type sequence._ |
| struct | [**CollectAllGridsOnDim**](structconnectivity__details_1_1CollectAllGridsOnDim.md) &lt;class StartPatch, class Grid1D, class InterfaceTypeSeq&gt;<br>_A class which collects all grids along a given dimension in both directions._  |
| struct | [**CollectAllInterfacesOnDim**](structconnectivity__details_1_1CollectAllInterfacesOnDim.md) &lt;class StartPatch, class Grid1D, class InterfaceTypeSeq&gt;<br>_A class which collects all grids along a given dimension in both directions._  |
| struct | [**CollectGridsAlongDim**](structconnectivity__details_1_1CollectGridsAlongDim.md) &lt;class StartEdge, class InterfaceTypeSeq, insert\_pos, class FoundGrids, class MatchingEdge, grid\_already\_found&gt;<br>_A class which collects grids along a given dimension on a specified direction from a starting edge._  |
| struct | [**CollectGridsAlongDim&lt; StartEdge, InterfaceTypeSeq, insert\_pos, FoundGrids, MatchingEdge, false &gt;**](structconnectivity__details_1_1CollectGridsAlongDim_3_01StartEdge_00_01InterfaceTypeSeq_00_01insd7a4bdb826ecb568487bbd509c5f008b.md) &lt;class StartEdge, class InterfaceTypeSeq, insert\_pos, class FoundGrids, class MatchingEdge&gt;<br>_Specialisation of_ [_**CollectGridsAlongDim**_](structconnectivity__details_1_1CollectGridsAlongDim.md) _to iterate recursively over the grids on the dimension._ |
| struct | [**CollectGridsAlongDim&lt; StartEdge, InterfaceTypeSeq, insert\_pos, FoundGrids, MatchingEdge, true &gt;**](structconnectivity__details_1_1CollectGridsAlongDim_3_01StartEdge_00_01InterfaceTypeSeq_00_01ins8ee738c554d8fbbf6bab92ba87dd3b80.md) &lt;class StartEdge, class InterfaceTypeSeq, insert\_pos, class FoundGrids, class MatchingEdge&gt;<br>_Specialisation of_ [_**CollectGridsAlongDim**_](structconnectivity__details_1_1CollectGridsAlongDim.md) _to stop when the grid has already been identified (due to periodicity)._ |
| struct | [**CollectGridsAlongDim&lt; StartEdge, InterfaceTypeSeq, insert\_pos, FoundGrids, OutsideEdge, false &gt;**](structconnectivity__details_1_1CollectGridsAlongDim_3_01StartEdge_00_01InterfaceTypeSeq_00_01ins70deef724c6e45ed62db534c3d9697ec.md) &lt;class StartEdge, class InterfaceTypeSeq, insert\_pos, class FoundGrids&gt;<br>_Specialisation of_ [_**CollectGridsAlongDim**_](structconnectivity__details_1_1CollectGridsAlongDim.md) _to stop when there are no more grids._ |
| struct | [**CollectInterfacesAlongDim**](structconnectivity__details_1_1CollectInterfacesAlongDim.md) &lt;class StartEdge, class InterfaceTypeSeq, insert\_pos, class FoundInterfaces, class MatchingEdge, interface\_already\_found&gt;<br>_A class which collects interfaces along a given dimension on a specified direction from a starting edge._  |
| struct | [**CollectInterfacesAlongDim&lt; StartEdge, InterfaceTypeSeq, insert\_pos, FoundInterfaces, MatchingEdge, false &gt;**](structconnectivity__details_1_1CollectInterfacesAlongDim_3_01StartEdge_00_01InterfaceTypeSeq_00_b2108f65f3430e895714f416a2f43701.md) &lt;class StartEdge, class InterfaceTypeSeq, insert\_pos, class FoundInterfaces, class MatchingEdge&gt;<br>_Specialisation of_ [_**CollectGridsAlongDim**_](structconnectivity__details_1_1CollectGridsAlongDim.md) _to iterate recursively over the grids on the dimension._ |
| struct | [**CollectInterfacesAlongDim&lt; StartEdge, InterfaceTypeSeq, insert\_pos, FoundInterfaces, MatchingEdge, true &gt;**](structconnectivity__details_1_1CollectInterfacesAlongDim_3_01StartEdge_00_01InterfaceTypeSeq_00_36879d7a164b5ac728612e4a981c6d65.md) &lt;class StartEdge, class InterfaceTypeSeq, insert\_pos, class FoundInterfaces, class MatchingEdge&gt;<br>_Specialisation of_ [_**CollectInterfacesAlongDim**_](structconnectivity__details_1_1CollectInterfacesAlongDim.md) _to stop when the interface has already been identified (due to periodicity)._ |
| struct | [**CollectInterfacesAlongDim&lt; StartEdge, InterfaceTypeSeq, insert\_pos, FoundInterfaces, OutsideEdge, false &gt;**](structconnectivity__details_1_1CollectInterfacesAlongDim_3_01StartEdge_00_01InterfaceTypeSeq_00_ebd86d7b2345baf351562d16964c47d9.md) &lt;class StartEdge, class InterfaceTypeSeq, insert\_pos, class FoundInterfaces&gt;<br>_Specialisation of_ [_**CollectGridsAlongDim**_](structconnectivity__details_1_1CollectGridsAlongDim.md) _to stop when there are no more grids._ |
| struct | [**EnforceFirstInterfaceEdge**](structconnectivity__details_1_1EnforceFirstInterfaceEdge.md) &lt;class InterfaceType, class FirstEdge&gt;<br>_A class to flip the edges in an interface to ensure that the correct edge comes first._  |
| struct | [**EnforceFirstInterfaceEdge&lt; Interface&lt; Edge2, FirstEdge, Orientations &gt;, FirstEdge &gt;**](structconnectivity__details_1_1EnforceFirstInterfaceEdge_3_01Interface_3_01Edge2_00_01FirstEdge_221a02b03250a49af1745b2263467420.md) &lt;class FirstEdge, class Edge2, Orientations&gt;<br>_Specialisation of_ [_**EnforceFirstInterfaceEdge**_](structconnectivity__details_1_1EnforceFirstInterfaceEdge.md) _for an interface which needs rearranging._ |
| struct | [**EnforceFirstInterfaceEdge&lt; Interface&lt; FirstEdge, Edge2, Orientations &gt;, FirstEdge &gt;**](structconnectivity__details_1_1EnforceFirstInterfaceEdge_3_01Interface_3_01FirstEdge_00_01Edge2_788676fcb3310ca4c1ec984ff0b4531b.md) &lt;class FirstEdge, class Edge2, Orientations&gt;<br>_Specialisation of_ [_**EnforceFirstInterfaceEdge**_](structconnectivity__details_1_1EnforceFirstInterfaceEdge.md) _for an interface which is already correctly arranged._ |
| struct | [**ExtractPatches**](structconnectivity__details_1_1ExtractPatches.md) &lt;class TypeSeq&gt;<br>_A class to find all the patches used by the various edges._  |
| struct | [**ExtractPatches&lt; ddc::detail::TypeSeq&lt; EdgeType1, EdgeTypes... &gt; &gt;**](structconnectivity__details_1_1ExtractPatches_3_01ddc_1_1detail_1_1TypeSeq_3_01EdgeType1_00_01EdgeTypes_8_8_8_01_4_01_4.md) &lt;class EdgeType1, EdgeTypes&gt;<br>_Specialisation of_ [_**ExtractPatches**_](structconnectivity__details_1_1ExtractPatches.md) _to iterate recursively over the edge type sequence._ |
| struct | [**ExtractPatches&lt; ddc::detail::TypeSeq&lt;&gt; &gt;**](structconnectivity__details_1_1ExtractPatches_3_01ddc_1_1detail_1_1TypeSeq_3_4_01_4.md) &lt;&gt;<br>_Specialisation of_ [_**ExtractPatches**_](structconnectivity__details_1_1ExtractPatches.md) _for an empty patch list._ |
| struct | [**FindInterface**](structconnectivity__details_1_1FindInterface.md) &lt;class [**Edge**](structEdge.md), class InterfaceTypeSeq&gt;<br>_A class to locate an interface which contains the specified edge._  |
| struct | [**FindInterface&lt; Edge, ddc::detail::TypeSeq&lt; Interface1, RemainingInterfaceTypes... &gt; &gt;**](structconnectivity__details_1_1FindInterface_3_01Edge_00_01ddc_1_1detail_1_1TypeSeq_3_01Interfac6d31b188ee73012ad6c98be99219379f.md) &lt;class [**Edge**](structEdge.md), class Interface1, RemainingInterfaceTypes&gt;<br>_Specialisation of_ [_**FindInterface**_](structconnectivity__details_1_1FindInterface.md) _to iterate recursively over the interface type sequence._ |
| struct | [**FindInterface&lt; Edge, ddc::detail::TypeSeq&lt; Interface&lt; Edge, OEdge, Orientations &gt;, RemainingInterfaceTypes... &gt; &gt;**](structconnectivity__details_1_1FindInterface_3_01Edge_00_01ddc_1_1detail_1_1TypeSeq_3_01Interfacd1aa547d7cc4bf022e85928246ab2d07.md) &lt;class [**Edge**](structEdge.md), class OEdge, Orientations, RemainingInterfaceTypes&gt;<br>_Specialisation of_ [_**FindInterface**_](structconnectivity__details_1_1FindInterface.md) _for the case where Edge1 from the first interface matches_[_**Edge**_](structEdge.md) _._ |
| struct | [**FindInterface&lt; Edge, ddc::detail::TypeSeq&lt; Interface&lt; OEdge, Edge, Orientations &gt;, RemainingInterfaceTypes... &gt; &gt;**](structconnectivity__details_1_1FindInterface_3_01Edge_00_01ddc_1_1detail_1_1TypeSeq_3_01Interfacee698732bdf35f06db097afe1714904c.md) &lt;class [**Edge**](structEdge.md), class OEdge, Orientations, RemainingInterfaceTypes&gt;<br>_Specialisation of_ [_**FindInterface**_](structconnectivity__details_1_1FindInterface.md) _for the case where Edge1 from the second interface matches_[_**Edge**_](structEdge.md) _._ |
| struct | [**FindInterface&lt; Edge, ddc::detail::TypeSeq&lt;&gt; &gt;**](structconnectivity__details_1_1FindInterface_3_01Edge_00_01ddc_1_1detail_1_1TypeSeq_3_4_01_4.md) &lt;class [**Edge**](structEdge.md)&gt;<br>_Specialisation of_ [_**FindInterface**_](structconnectivity__details_1_1FindInterface.md) _for an empty interface list._ |
| struct | [**FindPatch**](structconnectivity__details_1_1FindPatch.md) &lt;class Grid1D, class PatchTypeSeq&gt;<br>_A class to locate a patch which contains the specified grid._  |
| struct | [**FindPatch&lt; Grid1D, ddc::detail::TypeSeq&lt; Patch1, RemainingPatchTypes... &gt; &gt;**](structconnectivity__details_1_1FindPatch_3_01Grid1D_00_01ddc_1_1detail_1_1TypeSeq_3_01Patch1_00_33770856242f7c5cee1ce419b2efaf64.md) &lt;class Grid1D, class Patch1, RemainingPatchTypes&gt;<br>_Specialisation of_ [_**FindPatch**_](structconnectivity__details_1_1FindPatch.md) _to iterate recursively over the patch type sequence._ |
| struct | [**FindPatch&lt; Grid1D, ddc::detail::TypeSeq&lt;&gt; &gt;**](structconnectivity__details_1_1FindPatch_3_01Grid1D_00_01ddc_1_1detail_1_1TypeSeq_3_4_01_4.md) &lt;class Grid1D&gt;<br>_Specialisation of_ [_**FindPatch**_](structconnectivity__details_1_1FindPatch.md) _for an empty patch list._ |
| struct | [**FindPatch&lt; QueryGrid1D, ddc::detail::TypeSeq&lt; Patch&lt; OGrid, QueryGrid1D, BSpl1, BSpl2 &gt;, RemainingPatchTypes... &gt; &gt;**](structconnectivity__details_1_1FindPatch_3_01QueryGrid1D_00_01ddc_1_1detail_1_1TypeSeq_3_01Patch5f5acd76cfd59a22ebf513823679a320.md) &lt;class QueryGrid1D, class OGrid, class BSpl1, class BSpl2, RemainingPatchTypes&gt;<br> |
| struct | [**FindPatch&lt; QueryGrid1D, ddc::detail::TypeSeq&lt; Patch&lt; QueryGrid1D, OGrid, BSpl1, BSpl2 &gt;, RemainingPatchTypes... &gt; &gt;**](structconnectivity__details_1_1FindPatch_3_01QueryGrid1D_00_01ddc_1_1detail_1_1TypeSeq_3_01Patchd8fc8921dec760f8fe4c90c2a6947228.md) &lt;class QueryGrid1D, class OGrid, class BSpl1, class BSpl2, RemainingPatchTypes&gt;<br> |
| struct | [**FindRelevantIdxRangeType**](structconnectivity__details_1_1FindRelevantIdxRangeType.md) &lt;class QueryGrid1D, class IdxRangeTuple&gt;<br>_A class to find any index range types which contain an index range defined on the provided grid. E.g. Grid1, std::tuple&lt;IdxRange&lt;Grid1, Grid2&gt;, IdxRange&lt;Grid3,Grid4&gt;&gt; will find: ddc::detail::TypeSeq&lt;IdxRange&lt;Grid1, Grid2&gt;&gt;_  |
| struct | [**FindRelevantIdxRangeType&lt; QueryGrid1D, std::tuple&lt; IdxRangeHead, IdxRangeTypes... &gt; &gt;**](structconnectivity__details_1_1FindRelevantIdxRangeType_3_01QueryGrid1D_00_01std_1_1tuple_3_01Id3b131c802b30082f4412eb4689d6d53b.md) &lt;class QueryGrid1D, class IdxRangeHead, IdxRangeTypes&gt;<br>_Specialisation of_ [_**FindRelevantIdxRangeType**_](structconnectivity__details_1_1FindRelevantIdxRangeType.md) _to iterate recursively over the possible index range types._ |
| struct | [**FindRelevantIdxRangeType&lt; QueryGrid1D, std::tuple&lt;&gt; &gt;**](structconnectivity__details_1_1FindRelevantIdxRangeType_3_01QueryGrid1D_00_01std_1_1tuple_3_4_01_4.md) &lt;class QueryGrid1D&gt;<br>_Specialisation of_ [_**FindRelevantIdxRangeType**_](structconnectivity__details_1_1FindRelevantIdxRangeType.md) _for an empty list of index range types._ |
| struct | [**PatchConnection**](structconnectivity__details_1_1PatchConnection.md) &lt;class [**Patch**](structPatch.md), class InterfaceTypeSeq&gt;<br>_A class which finds all interfaces connected to a given patch._  |
| struct | [**PatchConnection&lt; Patch, ddc::detail::TypeSeq&lt; InterfaceType &gt; &gt;**](structconnectivity__details_1_1PatchConnection_3_01Patch_00_01ddc_1_1detail_1_1TypeSeq_3_01InterfaceType_01_4_01_4.md) &lt;class [**Patch**](structPatch.md), class InterfaceType&gt;<br>_Specialisation of_ [_**PatchConnection**_](structconnectivity__details_1_1PatchConnection.md) _for an interface list with one element._ |
| struct | [**PatchConnection&lt; Patch, ddc::detail::TypeSeq&lt; InterfaceType1, RemainingInterfaceTypes... &gt; &gt;**](structconnectivity__details_1_1PatchConnection_3_01Patch_00_01ddc_1_1detail_1_1TypeSeq_3_01Interd9a0a5e7aafe0b71fe7c76720b7c5da6.md) &lt;class [**Patch**](structPatch.md), class InterfaceType1, RemainingInterfaceTypes&gt;<br>_Specialisation of_ [_**PatchConnection**_](structconnectivity__details_1_1PatchConnection.md) _to iterate recursively over the interface type sequence._ |
| struct | [**PatchConnection&lt; Patch, ddc::detail::TypeSeq&lt;&gt; &gt;**](structconnectivity__details_1_1PatchConnection_3_01Patch_00_01ddc_1_1detail_1_1TypeSeq_3_4_01_4.md) &lt;class [**Patch**](structPatch.md)&gt;<br>_Specialisation of_ [_**PatchConnection**_](structconnectivity__details_1_1PatchConnection.md) _for an empty interface list._ |
| struct | [**SelectRelevantIdxRangeType**](structconnectivity__details_1_1SelectRelevantIdxRangeType.md) &lt;class QueryGrid1D, class IdxRangeType&gt;<br>_A class to create a type sequence which contains the index range if it can be used to index the grid._  |
| struct | [**SelectRelevantIdxRangeType&lt; QueryGrid1D, IdxRange&lt; IdxRangeGrids... &gt; &gt;**](structconnectivity__details_1_1SelectRelevantIdxRangeType_3_01QueryGrid1D_00_01IdxRange_3_01IdxRangeGrids_8_8_8_01_4_01_4.md) &lt;class QueryGrid1D, IdxRangeGrids&gt;<br>_Specialisation of_ [_**SelectRelevantIdxRangeType**_](structconnectivity__details_1_1SelectRelevantIdxRangeType.md) _to get access to the grids in the index range._ |
| struct | [**StripOutsideEdges**](structconnectivity__details_1_1StripOutsideEdges.md) &lt;class TypeSeq&gt;<br>_A class which finds all edges which are not_ [_**OutsideEdge**_](structOutsideEdge.md) _types._ |
| struct | [**StripOutsideEdges&lt; ddc::detail::TypeSeq&lt; EdgeType &gt; &gt;**](structconnectivity__details_1_1StripOutsideEdges_3_01ddc_1_1detail_1_1TypeSeq_3_01EdgeType_01_4_01_4.md) &lt;class EdgeType&gt;<br>_Specialisation of_ [_**StripOutsideEdges**_](structconnectivity__details_1_1StripOutsideEdges.md) _for the case with one edge in the list._ |
| struct | [**StripOutsideEdges&lt; ddc::detail::TypeSeq&lt; EdgeType1, RemainingEdgeTypes... &gt; &gt;**](structconnectivity__details_1_1StripOutsideEdges_3_01ddc_1_1detail_1_1TypeSeq_3_01EdgeType1_00_036e9ce7e4506982efa52c09ca049ae90.md) &lt;class EdgeType1, RemainingEdgeTypes&gt;<br>_Specialisation of_ [_**StripOutsideEdges**_](structconnectivity__details_1_1StripOutsideEdges.md) _to iterate recursively over the edge type sequence._ |
| struct | [**SwapExtremity**](structconnectivity__details_1_1SwapExtremity.md) &lt;class EdgeType&gt;<br>_A class to get the opposite edge of a grid line from one of the edges._  |
| struct | [**SwapExtremity&lt; Edge&lt; Patch, Grid1D, BACK &gt; &gt;**](structconnectivity__details_1_1SwapExtremity_3_01Edge_3_01Patch_00_01Grid1D_00_01BACK_01_4_01_4.md) &lt;class [**Patch**](structPatch.md), class Grid1D&gt;<br>_Specialisation of_ [_**SwapExtremity**_](structconnectivity__details_1_1SwapExtremity.md) _for an edge at the back end of a grid line._ |
| struct | [**SwapExtremity&lt; Edge&lt; Patch, Grid1D, FRONT &gt; &gt;**](structconnectivity__details_1_1SwapExtremity_3_01Edge_3_01Patch_00_01Grid1D_00_01FRONT_01_4_01_4.md) &lt;class [**Patch**](structPatch.md), class Grid1D&gt;<br>_Specialisation of_ [_**SwapExtremity**_](structconnectivity__details_1_1SwapExtremity.md) _for an edge at the front end of a grid line._ |
| struct | [**ToTuple**](structconnectivity__details_1_1ToTuple.md) &lt;class TypeSeq&gt;<br>_A class to convert a type sequence to a tuple type._  |
| struct | [**ToTuple&lt; ddc::detail::TypeSeq&lt; I... &gt; &gt;**](structconnectivity__details_1_1ToTuple_3_01ddc_1_1detail_1_1TypeSeq_3_01I_8_8_8_01_4_01_4.md) &lt;I&gt;<br>_Specialisation of_ [_**ToTuple**_](structconnectivity__details_1_1ToTuple.md) _for type sequences._ |


## Public Types

| Type | Name |
| ---: | :--- |
| enum  | [**InsertPosition**](#enum-insertposition)  <br>_An enumerate to help when inserting elements into a type sequence._  |
| typedef typename [**EnforceFirstInterfaceEdge**](structconnectivity__details_1_1EnforceFirstInterfaceEdge.md)&lt; InterfaceType, FirstEdge &gt;::type | [**enforce\_first\_interface\_edge\_t**](#typedef-enforce_first_interface_edge_t)  <br>_A tool to flip the edges in an interface to ensure that the correct edge comes first._  |
















































## Public Types Documentation




### enum InsertPosition 

_An enumerate to help when inserting elements into a type sequence._ 
```C++
enum connectivity_details::InsertPosition {
    FrontInsert,
    BackInsert
};
```




<hr>



### typedef enforce\_first\_interface\_edge\_t 

_A tool to flip the edges in an interface to ensure that the correct edge comes first._ 
```C++
using connectivity_details::enforce_first_interface_edge_t = typedef typename EnforceFirstInterfaceEdge<InterfaceType, FirstEdge>::type;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/connectivity/connectivity_details.hpp`

