// SPDX-License-Identifier: MIT

#include "constantrate.hpp"

ConstantRate::ConstantRate(double const value) : m_rate_value(value) {}

DFieldSpX ConstantRate::operator()(
        DFieldSpX rate,
        DConstFieldSpX density,
        DConstFieldSpX temperature) const
{
    ddc::parallel_fill(rate, m_rate_value);
    return rate;
}
