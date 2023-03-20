// SPDX-License-Identifier: MIT

#pragma once

#include <memory>

#include <geometry.hpp>

class IInterpolatorX
{
public:
    virtual ~IInterpolatorX() = default;

    virtual DSpanX operator()(DSpanX inout_data, ViewX<CoordX> coordinates) const = 0;
};

class InterpolatorXProxy : public IInterpolatorX
{
    std::unique_ptr<IInterpolatorX> m_impl;

public:
    InterpolatorXProxy(std::unique_ptr<IInterpolatorX>&& impl) : m_impl(std::move(impl)) {}

    ~InterpolatorXProxy() override = default;

    virtual DSpanX operator()(DSpanX const inout_data, ViewX<CoordX> const coordinates)
            const override
    {
        return (*m_impl)(inout_data, coordinates);
    }
};

class IPreallocatableInterpolatorX : public IInterpolatorX
{
public:
    ~IPreallocatableInterpolatorX() override = default;

    virtual InterpolatorXProxy preallocate() const = 0;
};
