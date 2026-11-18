

# Struct connectivity\_details::FindPatch&lt; QueryGrid1D, ddc::TypeSeq&lt; Patch&lt; PatchId, QueryGrid1D, OGrid, BSpl1, BSpl2 &gt;, RemainingPatchTypes... &gt; &gt;

**template &lt;class QueryGrid1D, std::size\_t PatchId, class OGrid, class BSpl1, class BSpl2, class... RemainingPatchTypes&gt;**



[**ClassList**](annotated.md) **>** [**connectivity\_details**](namespaceconnectivity__details.md) **>** [**FindPatch&lt; QueryGrid1D, ddc::TypeSeq&lt; Patch&lt; PatchId, QueryGrid1D, OGrid, BSpl1, BSpl2 &gt;, RemainingPatchTypes... &gt; &gt;**](structconnectivity__details_1_1FindPatch_3_01QueryGrid1D_00_01ddc_1_1TypeSeq_3_01Patch_3_01Patch25efc32586ee939744293d0c53ed9d67.md)



[More...](#detailed-description)

* `#include <connectivity_details.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef [**Patch**](structPatch.md)&lt; PatchId, QueryGrid1D, OGrid, BSpl1, BSpl2 &gt; | [**type**](#typedef-type)  <br>_The type found by the class._  |
















































## Detailed Description


Specialisation of [**FindPatch**](structconnectivity__details_1_1FindPatch.md) for the case where Grid1 from first patch in type sequence matches QueryGrid1D 


    
## Public Types Documentation




### typedef type 

_The type found by the class._ 
```C++
using connectivity_details::FindPatch< QueryGrid1D, ddc::TypeSeq< Patch< PatchId, QueryGrid1D, OGrid, BSpl1, BSpl2 >, RemainingPatchTypes... > >::type =  Patch<PatchId, QueryGrid1D, OGrid, BSpl1, BSpl2>;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/connectivity/connectivity_details.hpp`

