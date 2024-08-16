// SPDX-License-Identifier: MIT

#pragma once

#include "geometry.hpp"

/**
 * @brief An abstract class for initializing the equilibrium state of the distribution function.
 * The equilibrium state does not depend on spatial dimensions.
 */
class IEquilibrium
{
public:
    virtual ~IEquilibrium() = default;

    /**
     * @brief Operator for initializing a distribution function that does not depend on space.
     *
     * @param[in, out] allfequilibrium On input: the uninitialized distribution function.
     *                                 On output: the initialized distribution function.
     * @return The initialized distribution function.
     */

    virtual DFieldSpVxVy operator()(DFieldSpVxVy allfequilibrium) const = 0;
};
