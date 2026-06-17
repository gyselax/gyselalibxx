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
struct ElementwiseChoice;

} // namespace polar_foot_finder_details
