

# File utils\_patch\_locators.hpp

[**File List**](files.md) **>** [**connectivity**](dir_28b51abc9241105ab41b66c468e7d019.md) **>** [**utils\_patch\_locators.hpp**](utils__patch__locators_8hpp.md)

[Go to the documentation of this file](utils__patch__locators_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once

#include "onion_patch_locator.hpp"


template <class Tag>
struct is_onion_patch_locator : std::false_type
{
};

template <class T>
inline constexpr bool is_onion_patch_locator_v = is_onion_patch_locator<T>::value;


template <
        class MultipatchIdxRanges,
        class LogicalToPhysicalMapping,
        class PhysicalToLogicalMapping,
        class ExecSpace>
struct is_onion_patch_locator<OnionPatchLocator<
        MultipatchIdxRanges,
        LogicalToPhysicalMapping,
        PhysicalToLogicalMapping,
        ExecSpace>> : std::true_type
{
};
```


