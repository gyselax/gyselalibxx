// SPDX-License-Identifier: MIT

#pragma once

#include <geometry.hpp>

enum class RhsType { Source, Sink };

class IRightHandSide
{
public:
    virtual ~IRightHandSide() = default;

    virtual DSpanSpXVx operator()(DSpanSpXVx allfdistribu, double dt) const = 0;
};
