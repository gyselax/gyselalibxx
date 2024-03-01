// SPDX-License-Identifier: MIT

#include <memory>

#include <ddc/ddc.hpp>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "geometry.hpp"
#include "iadvectionvx.hpp"
#include "iadvectionx.hpp"
#include "species_info.hpp"
#include "splitvlasovsolver.hpp"

class MockAdvectionX : public IAdvectionSpatial<GeometryXVx, IDimX>
{
    using DDom = typename GeometryXVx::FdistribuDDom;

public:
    MockAdvectionX() = default;

    MOCK_METHOD(
            (device_t<ddc::ChunkSpan<double, DDom>>),
            CallOp,
            ((device_t<ddc::ChunkSpan<double, DDom>>)allfdistribu, double dt),
            (const));

    device_t<ddc::ChunkSpan<double, DDom>> operator()(
            device_t<ddc::ChunkSpan<double, DDom>> const allfdistribu,
            double const dt) const override
    {
        return this->CallOp(allfdistribu, dt);
    }
};

class MockAdvectionVx : public IAdvectionVelocity<GeometryXVx, IDimVx>
{
    using DDom = typename GeometryXVx::FdistribuDDom;
    using IDimX = typename GeometryXVx::SpatialDDom;

public:
    MockAdvectionVx() = default;

    MOCK_METHOD(
            (device_t<ddc::ChunkSpan<double, DDom>>),
            CallOp,
            // clang-format off
            ((device_t<ddc::ChunkSpan<double, DDom>>) allfdistribu,
             (device_t<ddc::ChunkSpan<const double, IDimX>>) efield,
             (double) dt),
            // clang-format on
            (const));

    device_t<ddc::ChunkSpan<double, DDom>> operator()(
            device_t<ddc::ChunkSpan<double, DDom>> const allfdistribu,
            device_t<ddc::ChunkSpan<const double, IDimX>> const efield,
            double const dt) const override
    {
        return this->CallOp(allfdistribu, efield, dt);
    }
};

using namespace ::testing;

TEST(SplitVlasovSolver, Ordering)
{
    IDomainSpXVx const dom(IndexSpXVx(0, 0, 0), IVectSpXVx(0, 0, 0));
    DFieldSpXVx fdistribu(dom);
    DSpanSpXVx const fdistribu_s(fdistribu);
    DFieldX const efield(ddc::select<IDimX>(dom));
    double const dt = 0.;

    MockAdvectionX const advec_x;
    MockAdvectionVx const advec_vx;
    SplitVlasovSolver const solver(advec_x, advec_vx);

    {
        InSequence s;

        EXPECT_CALL(advec_x, CallOp).WillOnce(Return(fdistribu_s));
        EXPECT_CALL(advec_vx, CallOp).WillOnce(Return(fdistribu_s));
        EXPECT_CALL(advec_x, CallOp).WillOnce(Return(fdistribu_s));
    }

    solver(fdistribu_s, efield, dt);
}
