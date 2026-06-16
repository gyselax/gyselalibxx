

# File elementwise\_choice.hpp

[**File List**](files.md) **>** [**advection**](dir_b90fde0f10c67a9aef841a6e6700f1f6.md) **>** [**polar\_foot\_finders**](dir_71ecba75f1675b0dc7d7b6463874ff90.md) **>** [**elementwise\_choice.hpp**](elementwise__choice_8hpp.md)

[Go to the documentation of this file](elementwise__choice_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

enum class FootFindingSpace { LOGICAL, PSEUDO_PHYSICAL, PHYSICAL };

enum class AdvectionFieldSpace { PHYSICAL, LOGICAL };

namespace polar_foot_finder_details {

template <
        FootFindingSpace FFSpace,
        AdvectionFieldSpace AFSpace,
        class GridR,
        class GridTheta,
        class IdxRangeOperator,
        class RThetaAdvectionEvaluator,
        class AdvecCoefField,
        class TimeStepperBuilder,
        concepts::Mapping LogicalToPhysicalMapping>
class ElementwiseChoice;

} // namespace polar_foot_finder_details
```


