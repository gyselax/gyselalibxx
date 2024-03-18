// Copyright (C) The DDC development team, see COPYRIGHT.md file
//
// SPDX-License-Identifier: MIT

#include <utility>

#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include "relocatable_device_code_initialization.hpp"

std::pair<ddc::Coordinate<rdc::DimX>, ddc::Coordinate<rdc::DimX>> read_from_device()
{
    rdc::DDomX const dom_x(rdc::DElemX(0), rdc::DVectX(2));
    ddc::Chunk allocation_d(dom_x, ddc::DeviceAllocator<double>());
    ddc::ChunkSpan const array = allocation_d.span_view();
    ddc::parallel_for_each(
            dom_x.take_first(rdc::DVectX(1)),
            KOKKOS_LAMBDA(rdc::DElemX const ix) { array(ix) = ddc::origin<rdc::DDimX>(); });
    ddc::parallel_for_each(
            dom_x.take_last(rdc::DVectX(1)),
            KOKKOS_LAMBDA(rdc::DElemX const ix) { array(ix) = ddc::step<rdc::DDimX>(); });
    ddc::Chunk allocation_h(dom_x, ddc::HostAllocator<double>());
    ddc::parallel_deepcopy(allocation_h, allocation_d);
    return std::pair<
            ddc::Coordinate<rdc::DimX>,
            ddc::Coordinate<rdc::DimX>>(allocation_h(rdc::DElemX(0)), allocation_h(rdc::DElemX(1)));
}

TEST(RelocatableDeviceCode, ReadFromDevice)
{
    ddc::Coordinate<rdc::DimX> const origin(-1.);
    double const step = 0.5;
    rdc::initialize_ddimx(origin, step);
    auto const& [origin_from_device, step_from_device] = read_from_device();
    EXPECT_DOUBLE_EQ(origin, origin_from_device);
    EXPECT_DOUBLE_EQ(step, step_from_device);
}
