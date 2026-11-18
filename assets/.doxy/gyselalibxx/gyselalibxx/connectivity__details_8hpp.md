

# File connectivity\_details.hpp



[**FileList**](files.md) **>** [**connectivity**](dir_28b51abc9241105ab41b66c468e7d019.md) **>** [**connectivity\_details.hpp**](connectivity__details_8hpp.md)

[Go to the source code of this file](connectivity__details_8hpp_source.md)



* `#include <array>`
* `#include <ddc/ddc.hpp>`
* `#include "ddc_aliases.hpp"`
* `#include "edge.hpp"`
* `#include "interface.hpp"`
* `#include "patch.hpp"`
* `#include "type_seq_tools.hpp"`













## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**connectivity\_details**](namespaceconnectivity__details.md) <br> |


## Classes

| Type | Name |
| ---: | :--- |
| struct | [**AddToTypeSeq&lt; ToInsert, TypeSeq, BackInsert &gt;**](structconnectivity__details_1_1AddToTypeSeq_3_01ToInsert_00_01TypeSeq_00_01BackInsert_01_4.md) &lt;class ToInsert, TypeSeq&gt;<br>_Specialisation of_ [_**AddToTypeSeq**_](structconnectivity__details_1_1AddToTypeSeq.md) _to add an element at the back of the type sequence._ |
| struct | [**AddToTypeSeq&lt; ToInsert, TypeSeq, FrontInsert &gt;**](structconnectivity__details_1_1AddToTypeSeq_3_01ToInsert_00_01TypeSeq_00_01FrontInsert_01_4.md) &lt;class ToInsert, TypeSeq&gt;<br>_Specialisation of_ [_**AddToTypeSeq**_](structconnectivity__details_1_1AddToTypeSeq.md) _to add an element at the front of the type sequence._ |
| struct | [**CollectAllGridsOnDim**](structconnectivity__details_1_1CollectAllGridsOnDim.md) &lt;StartPatch, Grid1D, InterfaceTypeSeq&gt;<br>_A class which collects all grids along a given dimension in both directions._  |
| struct | [**CollectAllInterfacesOnDim**](structconnectivity__details_1_1CollectAllInterfacesOnDim.md) &lt;StartPatch, Grid1D, InterfaceTypeSeq&gt;<br>_A class which collects all grids along a given dimension in both directions._  |
| struct | [**CollectGridsAlongDim&lt; StartEdge, InterfaceTypeSeq, insert\_pos, FoundGrids, MatchingEdge, false &gt;**](structconnectivity__details_1_1CollectGridsAlongDim_3_01StartEdge_00_01InterfaceTypeSeq_00_01insd7a4bdb826ecb568487bbd509c5f008b.md) &lt;StartEdge, InterfaceTypeSeq, insert\_pos, FoundGrids, MatchingEdge&gt;<br>_Specialisation of_ [_**CollectGridsAlongDim**_](structconnectivity__details_1_1CollectGridsAlongDim.md) _to iterate recursively over the grids on the dimension._ |
| struct | [**CollectGridsAlongDim&lt; StartEdge, InterfaceTypeSeq, insert\_pos, FoundGrids, MatchingEdge, true &gt;**](structconnectivity__details_1_1CollectGridsAlongDim_3_01StartEdge_00_01InterfaceTypeSeq_00_01ins8ee738c554d8fbbf6bab92ba87dd3b80.md) &lt;StartEdge, InterfaceTypeSeq, insert\_pos, FoundGrids, MatchingEdge&gt;<br>_Specialisation of_ [_**CollectGridsAlongDim**_](structconnectivity__details_1_1CollectGridsAlongDim.md) _to stop when the grid has already been identified (due to periodicity)._ |
| struct | [**CollectGridsAlongDim&lt; StartEdge, InterfaceTypeSeq, insert\_pos, FoundGrids, OutsideEdge, false &gt;**](structconnectivity__details_1_1CollectGridsAlongDim_3_01StartEdge_00_01InterfaceTypeSeq_00_01ins70deef724c6e45ed62db534c3d9697ec.md) &lt;StartEdge, InterfaceTypeSeq, insert\_pos, FoundGrids&gt;<br>_Specialisation of_ [_**CollectGridsAlongDim**_](structconnectivity__details_1_1CollectGridsAlongDim.md) _to stop when there are no more grids._ |
| struct | [**CollectInterfacesAlongDim&lt; StartEdge, InterfaceTypeSeq, insert\_pos, FoundInterfaces, MatchingEdge, false &gt;**](structconnectivity__details_1_1CollectInterfacesAlongDim_3_01StartEdge_00_01InterfaceTypeSeq_00_b2108f65f3430e895714f416a2f43701.md) &lt;StartEdge, InterfaceTypeSeq, insert\_pos, class FoundInterfaces, MatchingEdge&gt;<br>_Specialisation of_ [_**CollectInterfacesAlongDim**_](structconnectivity__details_1_1CollectInterfacesAlongDim.md) _to iterate recursively over the interfaces on the dimension._ |
| struct | [**CollectInterfacesAlongDim&lt; StartEdge, InterfaceTypeSeq, insert\_pos, FoundInterfaces, MatchingEdge, true &gt;**](structconnectivity__details_1_1CollectInterfacesAlongDim_3_01StartEdge_00_01InterfaceTypeSeq_00_36879d7a164b5ac728612e4a981c6d65.md) &lt;StartEdge, InterfaceTypeSeq, insert\_pos, class FoundInterfaces, MatchingEdge&gt;<br>_Specialisation of_ [_**CollectInterfacesAlongDim**_](structconnectivity__details_1_1CollectInterfacesAlongDim.md) _to stop when the interface has already been identified (due to periodicity)._ |
| struct | [**CollectInterfacesAlongDim&lt; StartEdge, InterfaceTypeSeq, insert\_pos, FoundInterfaces, OutsideEdge, false &gt;**](structconnectivity__details_1_1CollectInterfacesAlongDim_3_01StartEdge_00_01InterfaceTypeSeq_00_ebd86d7b2345baf351562d16964c47d9.md) &lt;StartEdge, InterfaceTypeSeq, insert\_pos, class FoundInterfaces&gt;<br>_Specialisation of_ [_**CollectInterfacesAlongDim**_](structconnectivity__details_1_1CollectInterfacesAlongDim.md) _to stop when there are no more interfaces._ |
| struct | [**EnforceFirstInterfaceEdge&lt;::Interface&lt; Edge2, FirstEdge, Orientations &gt;, FirstEdge &gt;**](structconnectivity__details_1_1EnforceFirstInterfaceEdge_3_1_1Interface_3_01Edge2_00_01FirstEdged69ce672ec5713276e61b0743e3351f4.md) &lt;FirstEdge, Edge2, Orientations&gt;<br>_Specialisation of_ [_**EnforceFirstInterfaceEdge**_](structconnectivity__details_1_1EnforceFirstInterfaceEdge.md) _for an interface which needs rearranging._ |
| struct | [**EnforceFirstInterfaceEdge&lt;::Interface&lt; FirstEdge, Edge2, Orientations &gt;, FirstEdge &gt;**](structconnectivity__details_1_1EnforceFirstInterfaceEdge_3_1_1Interface_3_01FirstEdge_00_01Edge250b3559d62afd5e1fb475939efa70e42.md) &lt;FirstEdge, Edge2, Orientations&gt;<br>_Specialisation of_ [_**EnforceFirstInterfaceEdge**_](structconnectivity__details_1_1EnforceFirstInterfaceEdge.md) _for an interface which is already correctly arranged._ |
| struct | [**ExtractPatches&lt; ddc::TypeSeq&lt; EdgeType1, EdgeTypes... &gt; &gt;**](structconnectivity__details_1_1ExtractPatches_3_01ddc_1_1TypeSeq_3_01EdgeType1_00_01EdgeTypes_8_8_8_01_4_01_4.md) &lt;EdgeType1, EdgeTypes&gt;<br>_Specialisation of_ [_**ExtractPatches**_](structconnectivity__details_1_1ExtractPatches.md) _to iterate recursively over the edge type sequence._ |
| struct | [**ExtractPatches&lt; ddc::TypeSeq&lt;&gt; &gt;**](structconnectivity__details_1_1ExtractPatches_3_01ddc_1_1TypeSeq_3_4_01_4.md) &lt;&gt;<br>_Specialisation of_ [_**ExtractPatches**_](structconnectivity__details_1_1ExtractPatches.md) _for an empty patch list._ |
| struct | [**FindInterface&lt; Edge, ddc::TypeSeq&lt; Interface1, RemainingInterfaceTypes... &gt; &gt;**](structconnectivity__details_1_1FindInterface_3_01Edge_00_01ddc_1_1TypeSeq_3_01Interface1_00_01Re9a613ed41c91d250bbda6f232579e462.md) &lt;Edge, class Interface1, RemainingInterfaceTypes&gt;<br>_Specialisation of_ [_**FindInterface**_](structconnectivity__details_1_1FindInterface.md) _to iterate recursively over the interface type sequence._ |
| struct | [**FindInterface&lt; Edge, ddc::TypeSeq&lt;::Interface&lt; Edge, OEdge, Orientations &gt;, RemainingInterfaceTypes... &gt; &gt;**](structconnectivity__details_1_1FindInterface_3_01Edge_00_01ddc_1_1TypeSeq_3_1_1Interface_3_01Edg7548375c29d9adc73538674e82ae95ac.md) &lt;Edge, OEdge, Orientations, RemainingInterfaceTypes&gt;<br>_Specialisation of_ [_**FindInterface**_](structconnectivity__details_1_1FindInterface.md) _for the case where Edge1 from the first interface matches_[_**Edge**_](structEdge.md) _._ |
| struct | [**FindInterface&lt; Edge, ddc::TypeSeq&lt;::Interface&lt; OEdge, Edge, Orientations &gt;, RemainingInterfaceTypes... &gt; &gt;**](structconnectivity__details_1_1FindInterface_3_01Edge_00_01ddc_1_1TypeSeq_3_1_1Interface_3_01OEdcf579b47e460dc35ee239300e37577bd.md) &lt;Edge, OEdge, Orientations, RemainingInterfaceTypes&gt;<br>_Specialisation of_ [_**FindInterface**_](structconnectivity__details_1_1FindInterface.md) _for the case where Edge1 from the second interface matches_[_**Edge**_](structEdge.md) _._ |
| struct | [**FindInterface&lt; Edge, ddc::TypeSeq&lt;&gt; &gt;**](structconnectivity__details_1_1FindInterface_3_01Edge_00_01ddc_1_1TypeSeq_3_4_01_4.md) &lt;Edge&gt;<br>_Specialisation of_ [_**FindInterface**_](structconnectivity__details_1_1FindInterface.md) _for an empty interface list._ |
| struct | [**FindPatch&lt; Grid1D, ddc::TypeSeq&lt; Patch1, RemainingPatchTypes... &gt; &gt;**](structconnectivity__details_1_1FindPatch_3_01Grid1D_00_01ddc_1_1TypeSeq_3_01Patch1_00_01RemainingPatchTypes_8_8_8_01_4_01_4.md) &lt;Grid1D, Patch1, RemainingPatchTypes&gt;<br> |
| struct | [**FindPatch&lt; Grid1D, ddc::TypeSeq&lt;&gt; &gt;**](structconnectivity__details_1_1FindPatch_3_01Grid1D_00_01ddc_1_1TypeSeq_3_4_01_4.md) &lt;Grid1D&gt;<br>_Specialisation of_ [_**FindPatch**_](structconnectivity__details_1_1FindPatch.md) _for an empty patch list._ |
| struct | [**FindPatchImpl&lt; Grid1D, ddc::TypeSeq&lt; Patch1, RemainingPatchTypes... &gt;, false &gt;**](structconnectivity__details_1_1FindPatchImpl_3_01Grid1D_00_01ddc_1_1TypeSeq_3_01Patch1_00_01Rema719fedd1c49728e1b817d795613a178c.md) &lt;class Grid1D, class Patch1, RemainingPatchTypes&gt;<br>_Specialisation of_ [_**FindPatchImpl**_](structconnectivity__details_1_1FindPatchImpl.md) _to iterate recursively over the patch type sequence._ |
| struct | [**FindPatchImpl&lt; Grid1D, ddc::TypeSeq&lt; Patch1, RemainingPatchTypes... &gt;, true &gt;**](structconnectivity__details_1_1FindPatchImpl_3_01Grid1D_00_01ddc_1_1TypeSeq_3_01Patch1_00_01Remaf6552449b60d10fbdc68783fc4df2a48.md) &lt;class Grid1D, class Patch1, RemainingPatchTypes&gt;<br>_Specialisation of_ [_**FindPatchImpl**_](structconnectivity__details_1_1FindPatchImpl.md) _for the case where the first patch contains the grid._ |
| struct | [**FindRelevantIdxRangeType&lt; QueryGrid1D, std::tuple&lt; IdxRangeHead, IdxRangeTypes... &gt; &gt;**](structconnectivity__details_1_1FindRelevantIdxRangeType_3_01QueryGrid1D_00_01std_1_1tuple_3_01Id3b131c802b30082f4412eb4689d6d53b.md) &lt;QueryGrid1D, IdxRangeHead, IdxRangeTypes&gt;<br>_Specialisation of_ [_**FindRelevantIdxRangeType**_](structconnectivity__details_1_1FindRelevantIdxRangeType.md) _to iterate recursively over the possible index range types._ |
| struct | [**FindRelevantIdxRangeType&lt; QueryGrid1D, std::tuple&lt;&gt; &gt;**](structconnectivity__details_1_1FindRelevantIdxRangeType_3_01QueryGrid1D_00_01std_1_1tuple_3_4_01_4.md) &lt;QueryGrid1D&gt;<br>_Specialisation of_ [_**FindRelevantIdxRangeType**_](structconnectivity__details_1_1FindRelevantIdxRangeType.md) _for an empty list of index range types._ |
| struct | [**PatchConnection&lt; Patch, ddc::TypeSeq&lt; InterfaceType &gt; &gt;**](structconnectivity__details_1_1PatchConnection_3_01Patch_00_01ddc_1_1TypeSeq_3_01InterfaceType_01_4_01_4.md) &lt;Patch, class InterfaceType&gt;<br>_Specialisation of_ [_**PatchConnection**_](structconnectivity__details_1_1PatchConnection.md) _for an interface list with one element._ |
| struct | [**PatchConnection&lt; Patch, ddc::TypeSeq&lt; InterfaceType1, RemainingInterfaceTypes... &gt; &gt;**](structconnectivity__details_1_1PatchConnection_3_01Patch_00_01ddc_1_1TypeSeq_3_01InterfaceType1_06c0add68931ef7e3252d237eeb2c9f3.md) &lt;Patch, class InterfaceType1, RemainingInterfaceTypes&gt;<br>_Specialisation of_ [_**PatchConnection**_](structconnectivity__details_1_1PatchConnection.md) _to iterate recursively over the interface type sequence._ |
| struct | [**PatchConnection&lt; Patch, ddc::TypeSeq&lt;&gt; &gt;**](structconnectivity__details_1_1PatchConnection_3_01Patch_00_01ddc_1_1TypeSeq_3_4_01_4.md) &lt;Patch&gt;<br>_Specialisation of_ [_**PatchConnection**_](structconnectivity__details_1_1PatchConnection.md) _for an empty interface list._ |
| struct | [**PatchIdsAreUnique&lt; ddc::TypeSeq&lt; Patches... &gt; &gt;**](structconnectivity__details_1_1PatchIdsAreUnique_3_01ddc_1_1TypeSeq_3_01Patches_8_8_8_01_4_01_4.md) &lt;Patches&gt;<br>_Specialisation of_ [_**PatchIdsAreUnique**_](structconnectivity__details_1_1PatchIdsAreUnique.md) _to unpack the patches._ |
| struct | [**SelectRelevantIdxRangeType&lt; QueryGrid1D, IdxRange&lt; IdxRangeGrids... &gt; &gt;**](structconnectivity__details_1_1SelectRelevantIdxRangeType_3_01QueryGrid1D_00_01IdxRange_3_01IdxRangeGrids_8_8_8_01_4_01_4.md) &lt;QueryGrid1D, IdxRangeGrids&gt;<br>_Specialisation of_ [_**SelectRelevantIdxRangeType**_](structconnectivity__details_1_1SelectRelevantIdxRangeType.md) _to get access to the grids in the index range._ |
| struct | [**StripOutsideEdges&lt; ddc::TypeSeq&lt; EdgeType &gt; &gt;**](structconnectivity__details_1_1StripOutsideEdges_3_01ddc_1_1TypeSeq_3_01EdgeType_01_4_01_4.md) &lt;EdgeType&gt;<br>_Specialisation of_ [_**StripOutsideEdges**_](structconnectivity__details_1_1StripOutsideEdges.md) _for the case with one edge in the list._ |
| struct | [**StripOutsideEdges&lt; ddc::TypeSeq&lt; EdgeType1, RemainingEdgeTypes... &gt; &gt;**](structconnectivity__details_1_1StripOutsideEdges_3_01ddc_1_1TypeSeq_3_01EdgeType1_00_01RemainingEdgeTypes_8_8_8_01_4_01_4.md) &lt;EdgeType1, RemainingEdgeTypes&gt;<br>_Specialisation of_ [_**StripOutsideEdges**_](structconnectivity__details_1_1StripOutsideEdges.md) _to iterate recursively over the edge type sequence._ |
| struct | [**SwapExtremity&lt; Edge&lt; Patch, Grid1D, BACK &gt; &gt;**](structconnectivity__details_1_1SwapExtremity_3_01Edge_3_01Patch_00_01Grid1D_00_01BACK_01_4_01_4.md) &lt;Patch, Grid1D&gt;<br>_Specialisation of_ [_**SwapExtremity**_](structconnectivity__details_1_1SwapExtremity.md) _for an edge at the back end of a grid line._ |
| struct | [**SwapExtremity&lt; Edge&lt; Patch, Grid1D, FRONT &gt; &gt;**](structconnectivity__details_1_1SwapExtremity_3_01Edge_3_01Patch_00_01Grid1D_00_01FRONT_01_4_01_4.md) &lt;Patch, Grid1D&gt;<br>_Specialisation of_ [_**SwapExtremity**_](structconnectivity__details_1_1SwapExtremity.md) _for an edge at the front end of a grid line._ |
| struct | [**ToTuple&lt; ddc::TypeSeq&lt; I... &gt; &gt;**](structconnectivity__details_1_1ToTuple_3_01ddc_1_1TypeSeq_3_01I_8_8_8_01_4_01_4.md) &lt;I&gt;<br>_Specialisation of_ [_**ToTuple**_](structconnectivity__details_1_1ToTuple.md) _for type sequences._ |


## Public Types

| Type | Name |
| ---: | :--- |
| typedef typename connectivity\_details:: CollectAllGridsOnDim&lt; StartPatch, Grid1D, InterfaceTypeSeq &gt;::type | [**collect\_grids\_on\_dim\_t**](#typedef-collect_grids_on_dim_t)  <br>_A tool to collect all grids along a given line including the specified grid._  |
| typedef typename connectivity\_details:: CollectAllInterfacesOnDim&lt; StartPatch, Grid1D, InterfaceTypeSeq &gt;::type | [**collect\_interfaces\_on\_dim\_t**](#typedef-collect_interfaces_on_dim_t)  <br>_A tool to collect all interfaces along a given line including the specified grid._  |
| typedef typename connectivity\_details:: FindInterface&lt; StartEdge, InterfaceTypeSeq &gt;::type::template OtherEdge&lt; StartEdge &gt; | [**equivalent\_edge\_t**](#typedef-equivalent_edge_t)  <br>_A utility to get the other edge in an interface._  |
| typedef typename [**connectivity\_details::ExtractPatches**](structconnectivity__details_1_1ExtractPatches.md)&lt; EdgeTypeSeq &gt;::type | [**extract\_patches\_t**](#typedef-extract_patches_t)  <br>_A tool to find all the patches used by the various edges._  |
| typedef typename [**connectivity\_details::FindInterface**](structconnectivity__details_1_1FindInterface.md)&lt; EdgeType, InterfaceTypeSeq &gt;::type | [**find\_associated\_interface\_t**](#typedef-find_associated_interface_t)  <br>_A tool to find the interface which contains a specified edge._  |
| typedef typename [**connectivity\_details::FindPatch**](structconnectivity__details_1_1FindPatch.md)&lt; Grid1D, PatchTypeSeq &gt;::type | [**find\_patch\_t**](#typedef-find_patch_t)  <br>_A tool to find a patch which contains the specified grid._  |
| typedef ddc::type\_seq\_element\_t&lt; 0, typename [**connectivity\_details::FindRelevantIdxRangeType**](structconnectivity__details_1_1FindRelevantIdxRangeType.md)&lt; QueryGrid1D, IdxRangeTuple &gt;::type &gt; | [**find\_relevant\_idx\_range\_t**](#typedef-find_relevant_idx_range_t)  <br>_A tool to find the first multi-D index range which contains a specific grid._  |
| typedef typename [**connectivity\_details::PatchConnection**](structconnectivity__details_1_1PatchConnection.md)&lt; StartPatch, InterfaceTypeSeq &gt;::type | [**interfaces\_of\_patch\_t**](#typedef-interfaces_of_patch_t)  <br>_A tool to find all interfaces directly connected to the start patch._  |
| typedef typename [**connectivity\_details::StripOutsideEdges**](structconnectivity__details_1_1StripOutsideEdges.md)&lt; ddc::TypeSeq&lt; EdgeType... &gt; &gt;::type | [**strip\_outside\_edges\_t**](#typedef-strip_outside_edges_t)  <br>_A tool to find all edges which are not_ [_**OutsideEdge**_](structOutsideEdge.md) _types._ |
| typedef typename [**connectivity\_details::ToTuple**](structconnectivity__details_1_1ToTuple.md)&lt; TypeSeq &gt;::type | [**to\_tuple\_t**](#typedef-to_tuple_t)  <br>_A tool for converting a DDC type sequence to a tuple type._  |




## Public Attributes

| Type | Name |
| ---: | :--- |
|  bool | [**patch\_ids\_are\_unique\_v**](#variable-patch_ids_are_unique_v)   = `[**connectivity\_details::PatchIdsAreUnique**](structconnectivity__details_1_1PatchIdsAreUnique.md)&lt;PatchTypeSeq&gt;::value`<br>_A tool to check that all the patches in a type sequence have different ids._  |












































## Public Types Documentation




### typedef collect\_grids\_on\_dim\_t 

_A tool to collect all grids along a given line including the specified grid._ 
```C++
using collect_grids_on_dim_t =  typename connectivity_details::
        CollectAllGridsOnDim<StartPatch, Grid1D, InterfaceTypeSeq>::type;
```




<hr>



### typedef collect\_interfaces\_on\_dim\_t 

_A tool to collect all interfaces along a given line including the specified grid._ 
```C++
using collect_interfaces_on_dim_t =  typename connectivity_details::
        CollectAllInterfacesOnDim<StartPatch, Grid1D, InterfaceTypeSeq>::type;
```




<hr>



### typedef equivalent\_edge\_t 

_A utility to get the other edge in an interface._ 
```C++
using equivalent_edge_t =  typename connectivity_details::
        FindInterface<StartEdge, InterfaceTypeSeq>::type::template OtherEdge<StartEdge>;
```





**Template parameters:**


* `StartEdge` The edge whose equivalent we are looking for. 
* `InterfaceTypeSeq` The DDC type sequence describing all possible interfaces which might contain the start edge. 




        

<hr>



### typedef extract\_patches\_t 

_A tool to find all the patches used by the various edges._ 
```C++
using extract_patches_t =  typename connectivity_details::ExtractPatches<EdgeTypeSeq>::type;
```




<hr>



### typedef find\_associated\_interface\_t 

_A tool to find the interface which contains a specified edge._ 
```C++
using find_associated_interface_t = 
        typename connectivity_details::FindInterface<EdgeType, InterfaceTypeSeq>::type;
```




<hr>



### typedef find\_patch\_t 

_A tool to find a patch which contains the specified grid._ 
```C++
using find_patch_t =  typename connectivity_details::FindPatch<Grid1D, PatchTypeSeq>::type;
```




<hr>



### typedef find\_relevant\_idx\_range\_t 

_A tool to find the first multi-D index range which contains a specific grid._ 
```C++
using find_relevant_idx_range_t =  ddc::type_seq_element_t<
        0,
        typename connectivity_details::FindRelevantIdxRangeType<QueryGrid1D, IdxRangeTuple>::type>;
```




<hr>



### typedef interfaces\_of\_patch\_t 

_A tool to find all interfaces directly connected to the start patch._ 
```C++
using interfaces_of_patch_t = 
        typename connectivity_details::PatchConnection<StartPatch, InterfaceTypeSeq>::type;
```




<hr>



### typedef strip\_outside\_edges\_t 

_A tool to find all edges which are not_ [_**OutsideEdge**_](structOutsideEdge.md) _types._
```C++
using strip_outside_edges_t = 
        typename connectivity_details::StripOutsideEdges<ddc::TypeSeq<EdgeType...>>::type;
```




<hr>



### typedef to\_tuple\_t 

_A tool for converting a DDC type sequence to a tuple type._ 
```C++
using to_tuple_t =  typename connectivity_details::ToTuple<TypeSeq>::type;
```




<hr>
## Public Attributes Documentation




### variable patch\_ids\_are\_unique\_v 

_A tool to check that all the patches in a type sequence have different ids._ 
```C++
bool patch_ids_are_unique_v;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/connectivity/connectivity_details.hpp`

