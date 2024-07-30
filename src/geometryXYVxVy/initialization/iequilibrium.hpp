// SPDX-License-Identifier: MIT

#pragma once

#include <geometry.hpp>

class IEquilibrium
{
public:
    virtual ~IEquilibrium() = default;

    virtual DSpanSpVxVy operator()(DSpanSpVxVy allfequilibrium) const = 0;
};
