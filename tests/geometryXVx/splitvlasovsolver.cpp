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

class MockAdvectionX : public IAdvectionSpatial<GeometryXVx, GridX>
{
    using IdxRange = typename GeometryXVx::IdxRangeFdistribu;

public:
    MockAdvectionX() = default;

    MOCK_METHOD((DField<IdxRange>), CallOp, ((DField<IdxRange>)allfdistribu, double dt), (const));

    DField<IdxRange> operator()(DField<IdxRange> const allfdistribu, double const dt) const override
    {
        return this->CallOp(allfdistribu, dt);
    }
};

class MockAdvectionVx : public IAdvectionVelocity<GeometryXVx, GridVx>
{
    using IdxRange = typename GeometryXVx::IdxRangeFdistribu;
    using GridX = typename GeometryXVx::IdxRangeSpatial;

public:
    MockAdvectionVx() = default;

    MOCK_METHOD(
            (DField<IdxRange>),
            CallOp,
            // clang-format off
            ((DField<IdxRange>) allfdistribu,
             (DConstField<GridX>) efield,
             (double) dt),
            // clang-format on
            (const));

    DField<IdxRange> operator()(
            DField<IdxRange> const allfdistribu,
            DConstField<GridX> const efield,
            double const dt) const override
    {
        return this->CallOp(allfdistribu, efield, dt);
    }
};

using namespace ::testing;

TEST(SplitVlasovSolver, Ordering)
{
    IdxRangeSpXVx const idx_range(IdxSpXVx(0, 0, 0), IdxStepSpXVx(0, 0, 0));
    DFieldMemSpXVx fdistribu(idx_range);
    DFieldSpXVx const fdistribu_s(fdistribu);
    DFieldMemX const efield(ddc::select<GridX>(idx_range));
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
