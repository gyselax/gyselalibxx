// SPDX-License-Identifier: MIT

#pragma once

#include <geometry.hpp>

class IEquilibrium
{
public:
    virtual ~IEquilibrium() = default;

    virtual device_t<DSpanSpVx> operator()(device_t<DSpanSpVx> allfequilibrium) const = 0;
};
