

# Struct connectivity\_details::FindPatches&lt; ddc::TypeSeq&lt; Grids... &gt;, PatchTypeSeq &gt;

**template &lt;ddc::concepts::discrete\_dimension... Grids, ddc::concepts::type\_seq PatchTypeSeq&gt;**



[**ClassList**](annotated.md) **>** [**connectivity\_details**](namespaceconnectivity__details.md) **>** [**FindPatches&lt; ddc::TypeSeq&lt; Grids... &gt;, PatchTypeSeq &gt;**](structconnectivity__details_1_1FindPatches_3_01ddc_1_1TypeSeq_3_01Grids_8_8_8_01_4_00_01PatchTypeSeq_01_4.md)



_Specialisation of_ [_**FindPatches**_](structconnectivity__details_1_1FindPatches.md) _to unpack the grids._

* `#include <connectivity_details.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef ddc::TypeSeq&lt; typename [**FindPatch**](structconnectivity__details_1_1FindPatch.md)&lt; Grids, PatchTypeSeq &gt;::type... &gt; | [**type**](#typedef-type)  <br>_The type found by the class._  |
















































## Public Types Documentation




### typedef type 

_The type found by the class._ 
```C++
using connectivity_details::FindPatches< ddc::TypeSeq< Grids... >, PatchTypeSeq >::type =  ddc::TypeSeq<typename FindPatch<Grids, PatchTypeSeq>::type...>;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/connectivity/connectivity_details.hpp`

