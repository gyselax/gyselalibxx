

# Struct connectivity\_details::PatchConnection&lt; Patch, ddc::detail::TypeSeq&lt; InterfaceType1, RemainingInterfaceTypes... &gt; &gt;

**template &lt;class [**Patch**](structPatch.md), class InterfaceType1, class... RemainingInterfaceTypes&gt;**



[**ClassList**](annotated.md) **>** [**connectivity\_details**](namespaceconnectivity__details.md) **>** [**PatchConnection&lt; Patch, ddc::detail::TypeSeq&lt; InterfaceType1, RemainingInterfaceTypes... &gt; &gt;**](structconnectivity__details_1_1PatchConnection_3_01Patch_00_01ddc_1_1detail_1_1TypeSeq_3_01Interd9a0a5e7aafe0b71fe7c76720b7c5da6.md)



_Specialisation of_ [_**PatchConnection**_](structconnectivity__details_1_1PatchConnection.md) _to iterate recursively over the interface type sequence._

* `#include <connectivity_details.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef ddc::type\_seq\_merge\_t&lt; typename [**PatchConnection**](structconnectivity__details_1_1PatchConnection.md)&lt; [**Patch**](structPatch.md), ddc::detail::TypeSeq&lt; InterfaceType1 &gt; &gt;::type, typename [**PatchConnection**](structconnectivity__details_1_1PatchConnection.md)&lt; [**Patch**](structPatch.md), ddc::detail::TypeSeq&lt; RemainingInterfaceTypes... &gt; &gt;::type &gt; | [**type**](#typedef-type)  <br>_The type found by the class._  |
















































## Public Types Documentation




### typedef type 

_The type found by the class._ 
```C++
using connectivity_details::PatchConnection< Patch, ddc::detail::TypeSeq< InterfaceType1, RemainingInterfaceTypes... > >::type =  ddc::type_seq_merge_t< typename PatchConnection<Patch, ddc::detail::TypeSeq<InterfaceType1> >::type, typename PatchConnection<Patch, ddc::detail::TypeSeq<RemainingInterfaceTypes...> >:: type>;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/connectivity/connectivity_details.hpp`

