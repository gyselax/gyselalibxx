#include <memory>

#include <ddc/BlockSpan>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "geometry.hpp"
#include "iadvectionvx.hpp"
#include "iadvectionx.hpp"
#include "species_info.hpp"
#include "splitvlasovsolver.hpp"

class MockAdvectionX : public IAdvectionX
{
public:
    MOCK_METHOD(DSpanSpXVx, CallOp, (DSpanSpXVx fdistribu, double dt), (const));
    DSpanSpXVx operator()(DSpanSpXVx fdistribu, double dt) const override
    {
        return this->CallOp(fdistribu, dt);
    }
};

class MockAdvectionVx : public IAdvectionVx
{
public:
    MOCK_METHOD(DSpanSpXVx, CallOp, (DSpanSpXVx fdistribu, DViewX efield, double dt), (const));
    DSpanSpXVx operator()(DSpanSpXVx fdistribu, DViewX efield, double dt) const override
    {
        return this->CallOp(fdistribu, efield, dt);
    }
};

using namespace ::testing;

TEST(SplitVlasovSolver, ordering)
{
    MeshX mesh_x(RCoordX(0.), RCoordX(2.));
    MeshVx mesh_vx(RCoordVx(0.), RCoordVx(2.));
    MeshSp mesh_sp;
    MDomainSpXVx const dom(mesh_sp, mesh_x, mesh_vx, MCoordSpXVx(0, 0, 0), MLengthSpXVx(0, 0, 0));
    DBlockSpXVx const fdistribu(dom);
    DSpanSpXVx const fdistribu_s(fdistribu);
    DBlockX const efield(select<MeshX>(dom));
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
