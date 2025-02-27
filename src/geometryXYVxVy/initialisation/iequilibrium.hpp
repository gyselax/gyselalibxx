// SPDX-License-Identifier: MIT

#pragma once

#include "geometry.hpp"

/**
 * @brief An abstract class for initialising the equilibrium state of the distribution function.
 * The equilibrium state does not depend on spatial dimensions.
 */
class IEquilibrium
{
public:
    virtual ~IEquilibrium() = default;

    /**
     * @brief Operator for initialising a distribution function that does not depend on space.
     *
     * @param[in, out] allfequilibrium On input: the uninitialized distribution function.
     *                                 On output: the initialised distribution function.
     * @return The initialised distribution function.
     */

    virtual DFieldSpVxVy operator()(DFieldSpVxVy allfequilibrium) const = 0;
};
