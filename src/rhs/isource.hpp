// SPDX-License-Identifier: MIT

#pragma once

#include <geometry.hpp>

class ISource
{
public:
    virtual ~ISource() = default;

    virtual DSpanSpXVx operator()(DSpanSpXVx allfdistribu, double dt) const = 0;
};
