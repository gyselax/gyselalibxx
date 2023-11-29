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
public:
    MockAdvectionX() = default;

    MOCK_METHOD(DSpanSpXVx, CallOp, (DSpanSpXVx fdistribu, double dt), (const));
    DSpanSpXVx operator()(DSpanSpXVx fdistribu, double dt) const override
    {
        return this->CallOp(fdistribu, dt);
    }
};

class MockAdvectionVx : public IAdvectionVelocity<GeometryXVx, IDimVx>
{
public:
    MockAdvectionVx() = default;

    MOCK_METHOD(DSpanSpXVx, CallOp, (DSpanSpXVx fdistribu, DViewX efield, double dt), (const));
    DSpanSpXVx operator()(DSpanSpXVx fdistribu, DViewX efield, double dt) const override
    {
        return this->CallOp(fdistribu, efield, dt);
    }
};

using namespace ::testing;

TEST(SplitVlasovSolver, Ordering)
{
    IDomainSpXVx const dom(IndexSpXVx(0, 0, 0), IVectSpXVx(0, 0, 0));
    device_t<DFieldSpXVx> fdistribu(dom);
    device_t<DSpanSpXVx> const fdistribu_s(fdistribu);
    DFieldX const efield(ddc::select<IDimX>(dom));
    double const dt = 0.;

    MockAdvectionX const advec_x;
    MockAdvectionVx const advec_vx;
    SplitVlasovSolver const solver(advec_x, advec_vx);

    auto fdistribu_s_host_alloc = ddc::create_mirror_view_and_copy(fdistribu_s);
    ddc::ChunkSpan fdistribu_s_host = fdistribu_s_host_alloc.span_view();
    {
        InSequence s;

        EXPECT_CALL(advec_x, CallOp).WillOnce(Return(fdistribu_s_host));
        EXPECT_CALL(advec_vx, CallOp).WillOnce(Return(fdistribu_s_host));
        EXPECT_CALL(advec_x, CallOp).WillOnce(Return(fdistribu_s_host));
    }

    solver(fdistribu_s, efield, dt);
}
