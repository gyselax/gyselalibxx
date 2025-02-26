// SPDX-License-Identifier: MIT

#pragma once

#include "geometry.hpp"

/**
 * @brief An abstract class that allows for initialising a distribution function.
 */
class IInitialisation
{
public:
    virtual ~IInitialisation() = default;

    /**
     * @brief Operator for initialising a distribution function.
     * @param[in, out] allfdistribu On input: the uninitialized distribution function.
     *                                 On output: the initialised distribution function.
     * @return The initialised distribution function.
     */
    virtual DFieldSpVparMu operator()(DFieldSpVparMu allfdistribu) const = 0;
};
