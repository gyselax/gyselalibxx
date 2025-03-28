

# File onion\_patch\_locator.hpp



[**FileList**](files.md) **>** [**connectivity**](dir_28b51abc9241105ab41b66c468e7d019.md) **>** [**onion\_patch\_locator.hpp**](onion__patch__locator_8hpp.md)

[Go to the source code of this file](onion__patch__locator_8hpp_source.md)



* `#include <stdexcept>`
* `#include <ddc/ddc.hpp>`
* `#include "ddc_aliases.hpp"`
* `#include "mapping_tools.hpp"`
* `#include "multipatch_type.hpp"`
* `#include "types.hpp"`















## Classes

| Type | Name |
| ---: | :--- |
| class | [**OnionPatchLocator&lt; MultipatchType&lt; IdxRangeOnPatch, Patches... &gt;, LogicalToPhysicalMapping, PhysicalToLogicalMapping, ExecSpace &gt;**](classOnionPatchLocator_3_01MultipatchType_3_01IdxRangeOnPatch_00_01Patches_8_8_8_01_4_00_01Logicff6c45b073183ccdfc0de0e4a415a7fa.md) &lt;Patches, class LogicalToPhysicalMapping, class PhysicalToLogicalMapping, class ExecSpace&gt;<br>[_**Patch**_](structPatch.md) _locator specialised for "onion" geometry._ |






















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**OnionPatchLocator**](#function-onionpatchlocator) (MultipatchIdxRanges const & all\_idx\_ranges, LogicalToPhysicalMapping const & to\_physical\_mapping, PhysicalToLogicalMapping const & to\_logical\_mapping) <br> |




























## Public Functions Documentation




### function OnionPatchLocator 

```C++
template<class MultipatchIdxRanges, class LogicalToPhysicalMapping, class PhysicalToLogicalMapping, class ExecSpace>
OnionPatchLocator (
    MultipatchIdxRanges const & all_idx_ranges,
    LogicalToPhysicalMapping const & to_physical_mapping,
    PhysicalToLogicalMapping const & to_logical_mapping
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/connectivity/onion_patch_locator.hpp`

