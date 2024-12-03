// SPDX-License-Identifier: MIT

#pragma once

#include "onion_patch_locator.hpp"


/// @brief Struct to identify if the patch locator is adapted to onion geometry.
template <class Tag>
struct is_onion_patch_locator : std::false_type
{
};

/// @brief Boolean: true if the patch locator is adapted to onion geometry.
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
