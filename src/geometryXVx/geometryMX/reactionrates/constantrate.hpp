// SPDX-License-Identifier: MIT

#pragma once

#include <geometry.hpp>

#include "ireactionrate.hpp"

/**
 * @brief A class that describes a reaction rate that is constant.
 */
class ConstantRate : public IReactionRate
{
private:
    double const m_rate_value;

public:
    /**
     * @brief Creates an instance of the ConstantRate class.
     * @param[in] value The constant value of the reaction rate. 
     */
    ConstantRate(double const value);

    ~ConstantRate() override = default;

    /**
     * @brief Operator for computing the reaction rate.
     * 
     * @param[out] rate On input: the uninitialized value of a Span referencing a reaction rate.
     *                      On output: the span referencing a constant reaction rate.
     * @param[in] density A span referencing the density.
     * @param[in] temperature A span referencing the temperatures.
     * @return A span referencing the reaction rate.
     */
    DSpanSpX operator()(DSpanSpX rate, DViewSpX density, DViewSpX temperature) const override;
};