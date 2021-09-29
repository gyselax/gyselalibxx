#include <memory>

#include <ddc/BlockSpan>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "geometry.h"
#include "iadvectionvx.h"
#include "iadvectionx.h"
#include "splitvlasovsolver.h"

class MockAdvectionX : public IAdvectionX
{
public:
    MOCK_METHOD(
            DSpanXVx,
            CallOp,
            (DSpanXVx fdistribu, double sqrt_me_on_mspecies, double dt),
            (const));
    DSpanXVx operator()(DSpanXVx fdistribu, double sqrt_me_on_mspecies, double dt) const override
    {
        return this->CallOp(fdistribu, sqrt_me_on_mspecies, dt);
    }
};

class MockAdvectionVx : public IAdvectionVx
{
public:
    MOCK_METHOD(
            DSpanXVx,
            CallOp,
            (DSpanXVx fdistribu, DViewX efield, double sqrt_me_on_mspecies, double dt),
            (const));
    DSpanXVx operator()(DSpanXVx fdistribu, DViewX efield, double sqrt_me_on_mspecies, double dt)
            const override
    {
        return this->CallOp(fdistribu, efield, sqrt_me_on_mspecies, dt);
    }
};

using namespace ::testing;

TEST(SplitVlasovSolver, ordering)
{
    MeshX mesh_x(RCoordX(0.), RCoordX(2.));
    MeshVx mesh_vx(RCoordVx(0.), RCoordVx(2.));
    ProductMesh mesh_x_vx(mesh_x, mesh_vx);
    MDomainXVx const dom(mesh_x_vx, MCoordXVx(0, 0), MCoordXVx(0, 0));
    DBlockXVx const fdistribu(dom);
    DSpanXVx const fdistribu_s(fdistribu);
    DBlockX const efield(ProductMDomain(get<MeshX>(dom)));
    double const sqrt_me_on_mspecies = 1.;
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

    solver(fdistribu_s, efield, sqrt_me_on_mspecies, dt);
}
