

# Struct connectivity\_details::FindPatch&lt; Grid1D, ddc::TypeSeq&lt; Patch1, RemainingPatchTypes... &gt; &gt;

**template &lt;class Grid1D, class Patch1, class... RemainingPatchTypes&gt;**



[**ClassList**](annotated.md) **>** [**connectivity\_details**](namespaceconnectivity__details.md) **>** [**FindPatch&lt; Grid1D, ddc::TypeSeq&lt; Patch1, RemainingPatchTypes... &gt; &gt;**](structconnectivity__details_1_1FindPatch_3_01Grid1D_00_01ddc_1_1TypeSeq_3_01Patch1_00_01RemainingPatchTypes_8_8_8_01_4_01_4.md)



_Specialisation of_ [_**FindPatch**_](structconnectivity__details_1_1FindPatch.md) _to iterate recursively over the patch type sequence._

* `#include <connectivity_details.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef typename [**FindPatch**](structconnectivity__details_1_1FindPatch.md)&lt; Grid1D, ddc::TypeSeq&lt; RemainingPatchTypes... &gt; &gt;::type | [**type**](#typedef-type)  <br>_The type found by the class._  |
















































## Public Types Documentation




### typedef type 

_The type found by the class._ 
```C++
using connectivity_details::FindPatch< Grid1D, ddc::TypeSeq< Patch1, RemainingPatchTypes... > >::type =  typename FindPatch<Grid1D, ddc::TypeSeq<RemainingPatchTypes...>>::type;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/connectivity/connectivity_details.hpp`

