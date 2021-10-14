#pragma once

#include <geometry.h>

class IInitialization
{
public:
    IInitialization() = default;

    IInitialization(IInitialization const& x) = default;

    IInitialization(IInitialization&& x) = default;

    virtual ~IInitialization() = default;

    IInitialization& operator=(IInitialization const& x) = default;

    IInitialization& operator=(IInitialization&& x) = default;

    virtual DSpanSpXVx operator()(DSpanSpXVx fdistribu) const = 0;
};
