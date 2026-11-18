

# Struct connectivity\_details::FindPatch&lt; Grid1D, ddc::TypeSeq&lt; Patch1, RemainingPatchTypes... &gt; &gt;

**template &lt;ddc::concepts::discrete\_dimension Grid1D, concepts::Patch Patch1, concepts::Patch... RemainingPatchTypes&gt;**



[**ClassList**](annotated.md) **>** [**connectivity\_details**](namespaceconnectivity__details.md) **>** [**FindPatch&lt; Grid1D, ddc::TypeSeq&lt; Patch1, RemainingPatchTypes... &gt; &gt;**](structconnectivity__details_1_1FindPatch_3_01Grid1D_00_01ddc_1_1TypeSeq_3_01Patch1_00_01RemainingPatchTypes_8_8_8_01_4_01_4.md)



[More...](#detailed-description)

* `#include <connectivity_details.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef typename [**FindPatchImpl**](structconnectivity__details_1_1FindPatchImpl.md)&lt; Grid1D, ddc::TypeSeq&lt; Patch1, RemainingPatchTypes... &gt;, ddc::in\_tags\_v&lt; Grid1D, ddc::TypeSeq&lt; typename Patch1::Grid1, typename Patch1::Grid2 &gt; &gt; &gt;:: [**type**](structconnectivity__details_1_1FindPatch_3_01Grid1D_00_01ddc_1_1TypeSeq_3_01Patch1_00_01RemainingPatchTypes_8_8_8_01_4_01_4.md#typedef-type) | [**type**](#typedef-type)  <br>_The type found by the class._  |
















































## Detailed Description


Specialisation of [**FindPatch**](structconnectivity__details_1_1FindPatch.md) which checks whether Grid1 or Grid2 of the first patch in the type sequence matches Grid1D, and otherwise iterates recursively over the patch type sequence. 


    
## Public Types Documentation




### typedef type 

_The type found by the class._ 
```C++
using connectivity_details::FindPatch< Grid1D, ddc::TypeSeq< Patch1, RemainingPatchTypes... > >::type =  typename FindPatchImpl<
            Grid1D,
            ddc::TypeSeq<Patch1, RemainingPatchTypes...>,
            ddc::in_tags_v<Grid1D, ddc::TypeSeq<typename Patch1::Grid1, typename Patch1::Grid2>>>::
            type;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/connectivity/connectivity_details.hpp`

