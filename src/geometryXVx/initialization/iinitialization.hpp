// SPDX-License-Identifier: MIT

#pragma once

#include <geometry.hpp>

class IInitialization
{
public:
    virtual ~IInitialization() = default;

    virtual device_t<DSpanSpXVx> operator()(device_t<DSpanSpXVx> allfdistribu) const = 0;
};
