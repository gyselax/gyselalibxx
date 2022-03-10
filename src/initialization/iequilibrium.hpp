// SPDX-License-Identifier: MIT

#pragma once

#include <geometry.hpp>

class IEquilibrium
{
public:
    virtual ~IEquilibrium() = default;

    virtual DSpanSpVx operator()(DSpanSpVx allfequilibrium) const = 0;
};
