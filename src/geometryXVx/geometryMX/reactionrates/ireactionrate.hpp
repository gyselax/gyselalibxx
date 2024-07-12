// SPDX-License-Identifier: MIT

#pragma once

#include <sll/view.hpp>

#include <geometry.hpp>

/**
 * @brief An abstract interface representing a reaction rate that depends on temperature and density.
 */
class IReactionRate
{
public:
    virtual ~IReactionRate() = default;

    /**
     * @brief Operator for computing the reaction rate from two Spans referencing density and temperature.
     * 
     * @param[out] rate On input: the uninitialized value of a Span referencing a reaction rate.
     *                      On output: the value of reaction rate.
     * @param[in] density The density at which the reaction rate should be computed.
     * @param[in] temperature The temperature at which the reaction rate should be computed.
     * @return A span referencing the reaction rate.
     */
    virtual DSpanSpX operator()(DSpanSpX rate, DViewSpX density, DViewSpX temperature) const = 0;
};
