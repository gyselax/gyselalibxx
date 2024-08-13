// SPDX-License-Identifier: MIT
#pragma once

#include <sll/view.hpp>

#include "geometry.hpp"

/**
 * @brief An abstract interface representing a reaction rate that depends on temperature and density.
 */
class IReactionRate
{
public:
    virtual ~IReactionRate() = default;

    /**
     * @brief Operator for computing the reaction rate from two Fields referencing density and temperature.
     * 
     * @param[out] rate On input: the uninitialized value of a Field referencing a reaction rate.
     *                      On output: the value of reaction rate.
     * @param[in] density The density at which the reaction rate should be computed.
     * @param[in] temperature The temperature at which the reaction rate should be computed.
     * @return A field referencing the reaction rate.
     */
    virtual DFieldSpX operator()(DFieldSpX rate, DConstFieldSpX density, DConstFieldSpX temperature)
            const = 0;
};
