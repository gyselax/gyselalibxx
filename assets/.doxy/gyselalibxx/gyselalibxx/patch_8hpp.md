

# File patch.hpp



[**FileList**](files.md) **>** [**connectivity**](dir_28b51abc9241105ab41b66c468e7d019.md) **>** [**patch.hpp**](patch_8hpp.md)

[Go to the source code of this file](patch_8hpp_source.md)



* `#include <type_traits>`
* `#include <utility>`
* `#include <ddc/ddc.hpp>`
* `#include "ddc_aliases.hpp"`













## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**concepts**](namespaceconcepts.md) <br> |


## Classes

| Type | Name |
| ---: | :--- |
| struct | [**Patch&lt; PatchId, grid1, grid2, bsplines\_dim1, bsplines\_dim2 &gt;**](structPatch_3_01PatchId_00_01grid1_00_01grid2_00_01bsplines__dim1_00_01bsplines__dim2_01_4.md) &lt;PatchId, class grid1, class grid2, class bsplines\_dim1, class bsplines\_dim2&gt;<br>_Tag for a patch._  |
| struct | [**is\_patch**](structis__patch.md) &lt;class PatchType&gt;<br>_A type trait to determine whether a type is a_ [_**Patch**_](structPatch.md) _, i.e. either an instantiation of_[_**Patch**_](structPatch.md) _or a class deriving from such an instantiation._ |






## Public Attributes

| Type | Name |
| ---: | :--- |
|  bool constexpr | [**is\_patch\_v**](#variable-is_patch_v)   = `[**is\_patch**](structis__patch.md)&lt;PatchType&gt;::value`<br> |












































## Public Attributes Documentation




### variable is\_patch\_v 

```C++
bool constexpr is_patch_v;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/connectivity/patch.hpp`

