// SPDX-License-Identifier: MIT

#pragma once

#include <memory>

#include <ddc/ddc.hpp>

template <class DDim>
class IInterpolator
{
    using CDim = typename DDim::continuous_dimension_type;

public:
    virtual ~IInterpolator() = default;

    virtual ddc::ChunkSpan<double, ddc::DiscreteDomain<DDim>> operator()(
            ddc::ChunkSpan<double, ddc::DiscreteDomain<DDim>> inout_data,
            ddc::ChunkSpan<const ddc::Coordinate<CDim>, ddc::DiscreteDomain<DDim>> coordinates)
            const = 0;
};

template <class DDim>
class InterpolatorProxy : public IInterpolator<DDim>
{
    using CDim = typename DDim::continuous_dimension_type;

    std::unique_ptr<IInterpolator<DDim>> m_impl;

public:
    InterpolatorProxy(std::unique_ptr<IInterpolator<DDim>>&& impl) : m_impl(std::move(impl)) {}

    ~InterpolatorProxy() override = default;

    ddc::ChunkSpan<double, ddc::DiscreteDomain<DDim>> operator()(
            ddc::ChunkSpan<double, ddc::DiscreteDomain<DDim>> inout_data,
            ddc::ChunkSpan<const ddc::Coordinate<CDim>, ddc::DiscreteDomain<DDim>> coordinates)
            const override
    {
        return (*m_impl)(inout_data, coordinates);
    }
};

template <class DDim>
class IPreallocatableInterpolator : public IInterpolator<DDim>
{
public:
    ~IPreallocatableInterpolator() override = default;

    virtual InterpolatorProxy<DDim> preallocate() const = 0;
};
