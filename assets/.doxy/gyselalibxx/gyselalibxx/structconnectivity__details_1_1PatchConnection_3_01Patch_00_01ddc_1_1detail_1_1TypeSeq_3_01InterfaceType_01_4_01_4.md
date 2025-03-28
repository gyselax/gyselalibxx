

# Struct connectivity\_details::PatchConnection&lt; Patch, ddc::detail::TypeSeq&lt; InterfaceType &gt; &gt;

**template &lt;class [**Patch**](structPatch.md), class InterfaceType&gt;**



[**ClassList**](annotated.md) **>** [**connectivity\_details**](namespaceconnectivity__details.md) **>** [**PatchConnection&lt; Patch, ddc::detail::TypeSeq&lt; InterfaceType &gt; &gt;**](structconnectivity__details_1_1PatchConnection_3_01Patch_00_01ddc_1_1detail_1_1TypeSeq_3_01InterfaceType_01_4_01_4.md)



_Specialisation of_ [_**PatchConnection**_](structconnectivity__details_1_1PatchConnection.md) _for an interface list with one element._

* `#include <connectivity_details.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef std::conditional\_t&lt; InterfaceType::template connected\_to\_patch&lt; [**Patch**](structPatch.md) &gt;(), ddc::detail::TypeSeq&lt; InterfaceType &gt;, ddc::detail::TypeSeq&lt;&gt; &gt; | [**type**](#typedef-type)  <br>_The type found by the class._  |
















































## Public Types Documentation




### typedef type 

_The type found by the class._ 
```C++
using connectivity_details::PatchConnection< Patch, ddc::detail::TypeSeq< InterfaceType > >::type =  std::conditional_t< InterfaceType::template connected_to_patch<Patch>(), ddc::detail::TypeSeq<InterfaceType>, ddc::detail::TypeSeq<> >;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/connectivity/connectivity_details.hpp`

