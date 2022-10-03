// SPDX-License-Identifier: MIT
#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include "relocatable_device_code_initialization.hpp"

std::pair<Coordinate<rdc::DimX>, Coordinate<rdc::DimX>> read_from_device()
{
    rdc::DDomX const dom_x(rdc::DElemX(0), rdc::DVectX(2));
    Chunk allocation_d(dom_x, DeviceAllocator<double>());
    ChunkSpan const array = allocation_d.span_view();
    for_each(
            policies::parallel_device,
            dom_x.take_first(rdc::DVectX(1)),
            DDC_LAMBDA(rdc::DElemX const ix) { array(ix) = origin<rdc::DDimX>(); });
    for_each(
            policies::parallel_device,
            dom_x.take_last(rdc::DVectX(1)),
            DDC_LAMBDA(rdc::DElemX const ix) { array(ix) = step<rdc::DDimX>(); });
    Chunk allocation_h(dom_x, HostAllocator<double>());
    deepcopy(allocation_h, allocation_d);
    return std::pair<
            Coordinate<rdc::DimX>,
            Coordinate<rdc::DimX>>(allocation_h(rdc::DElemX(0)), allocation_h(rdc::DElemX(1)));
}

TEST(RelocatableDeviceCode, ReadFromDevice)
{
    Coordinate<rdc::DimX> const origin(-1.);
    Coordinate<rdc::DimX> const step(0.5);
    rdc::initialize_ddimx(origin, step);
    auto const& [origin_from_device, step_from_device] = read_from_device();
    EXPECT_DOUBLE_EQ(origin, origin_from_device);
    EXPECT_DOUBLE_EQ(step, step_from_device);
}
