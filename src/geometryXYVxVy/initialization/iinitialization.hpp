// SPDX-License-Identifier: MIT

#pragma once

#include "geometry.hpp"

/**
 * @brief An abstract class that allows for initializing a distribution function.
 */
class IInitialization
{
public:
    virtual ~IInitialization() = default;

    /**
     * @brief Operator for initializing a distribution function.
     *
     * @param[in, out] allfdistribu On input: the uninitialized distribution function.
     *                                 On output: the initialized distribution function.
     * @return The initialized distribution function.
     */

    virtual DFieldSpXYVxVy operator()(DFieldSpXYVxVy allfdistribu) const = 0;
};
