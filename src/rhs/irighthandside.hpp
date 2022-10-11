// SPDX-License-Identifier: MIT

#pragma once

#include <geometry.hpp>

class IRightHandSide
{
public:
    virtual ~IRightHandSide() = default;

    virtual DSpanSpXVx operator()(DSpanSpXVx allfdistribu, double dt) const = 0;
};
