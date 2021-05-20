#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "block.h"
#include "blockview.h"
#include "splitvlasovsolver.h"

class MockAdvectionX : public IAdvectionX
{
public:
    MOCK_METHOD(
            DBlockSpanXVx,
            CallOp,
            (DBlockSpanXVx fdistribu, double mass_ratio, double dt),
            (const));
    DBlockSpanXVx operator()(DBlockSpanXVx fdistribu, double mass_ratio, double dt) const override
    {
        return this->CallOp(fdistribu, mass_ratio, dt);
    }
};

class MockAdvectionVx : public IAdvectionVx
{
public:
    MOCK_METHOD(
            DBlockSpanXVx,
            CallOp,
            (DBlockSpanXVx fdistribu, double mass_ratio, double dt),
            (const));
    DBlockSpanXVx operator()(DBlockSpanXVx fdistribu, double mass_ratio, double dt) const override
    {
        return this->CallOp(fdistribu, mass_ratio, dt);
    }
};

using namespace ::testing;

TEST(SplitVlasovSolver, ordering)
{
    MDomainXVx const dom(RCoordXVx(0., 0.), RCoordXVx(2., 2.), MCoordXVx(0, 0), MCoordXVx(1, 1));
    DBlockXVx const fdistribu(dom);
    DBlockSpanXVx const fdistribu_s(fdistribu);
    double const mass_ratio = 1.;
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

    solver(fdistribu_s, mass_ratio, dt);
}
