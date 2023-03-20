// SPDX-License-Identifier: MIT

#pragma once

#include <geometry.hpp>

class IInitialization
{
public:
    virtual ~IInitialization() = default;

    virtual DSpanSpXVx operator()(DSpanSpXVx allfdistribu) const = 0;
};
