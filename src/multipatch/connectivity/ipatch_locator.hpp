// SPDX-License-Identifier: MIT

#pragma once
#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"

/**
 * @brief Base class to locate on which patch a given physical coordinate is. 
 * 
 * The physical coordinates are supposed to be global 
 * and the logical coordinates are supposed to be locally defined.
 */
class IPatchLocator
{
public:
    /// @brief Default value to define outside domain (not a patch).
    static constexpr int outside_domain = -1;
};