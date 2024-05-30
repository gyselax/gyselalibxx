// SPDX-License-Identifier: MIT

#include "constantrate.hpp"

ConstantRate::ConstantRate(double const value) : m_rate_value(value) {}

DSpanSpX ConstantRate::operator()(DSpanSpX rate, DViewSpX density, DViewSpX temperature) const
{
    ddc::parallel_fill(rate, m_rate_value);
    return rate;
}
