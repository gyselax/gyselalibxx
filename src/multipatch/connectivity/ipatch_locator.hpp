// SPDX-License-Identifier: MIT

#pragma once
#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"

/**
 * @brief Base class to locate on which patch a given physical coordinate is. 
 * 
 * The physical coordinates are supposed to be global 
 * and the logical coordinates are supposed to be locally defined.
 * 
 * @tparam X First physical coordinate. 
 * @tparam Y Second physical coordinate. 
 */
template <class X, class Y>
class IPatchLocator
{
public:
    /// @brief Default value to define outside domain (not a patch).
    static constexpr int outside_domain = -1;

    IPatchLocator() = default;
    virtual ~IPatchLocator() = default;

    /**
     * @brief Get the patch index where the physical coordinate is. 
     * @param coord Physical coordinate we want to identify the associated patch. 
     * @return[int] Patch index of the associated patch. 
     */
    virtual int operator()(Coord<X, Y> const& coord) const = 0;
};