// SPDX-License-Identifier: MIT

#pragma once

#include <geometry.hpp>

class IInitialization
{
public:
    virtual ~IInitialization() = default;

    virtual DSpanSpXYVxVy operator()(DSpanSpXYVxVy allfdistribu) const = 0;
};
