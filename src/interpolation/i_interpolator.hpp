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
class IPreallocatableInterpolator : public IInterpolator<DDim>
{
public:
    ~IPreallocatableInterpolator() override = default;

    virtual std::unique_ptr<IInterpolator<DDim>> preallocate() const = 0;
};
