// SPDX-License-Identifier: MIT

#pragma once

#include "geometry.hpp"

/**
 * @brief An abstract class for initialising a distribution function that does not depend on space.
 */
class IEquilibrium
{
public:
    virtual ~IEquilibrium() = default;

    /**
     * @brief Operator for initialising a distribution function that does not depend on space.
     * @param[in, out] allfequilibrium On input: the uninitialized distribution function.
     *                                 On output: the initialised distribution function.
     * @return The initialised distribution function.
     */
    virtual DFieldSpVx operator()(DFieldSpVx allfequilibrium) const = 0;
};
