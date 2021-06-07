#include <gtest/gtest.h>

#include "block.h"

using namespace std;
using namespace std::experimental;

class DimX;
class DimVx;
using MeshX = UniformMesh<DimX>;
using MDomainX = UniformMDomain<DimX>;
using DBlockX = Block<MDomainX, double>;
using MeshXVx = UniformMesh<DimX, DimVx>;
using RCoordXVx = RCoord<DimX, DimVx>;
using MCoordXVx = MCoord<DimX, DimVx>;
using MDomainXVx = UniformMDomain<DimX, DimVx>;
using DBlockXVx = Block<MDomainXVx, double>;
using MeshVxX = UniformMesh<DimVx, DimX>;
using RCoordVxX = RCoord<DimVx, DimX>;
using MCoordVxX = MCoord<DimVx, DimX>;
using MDomainVxX = UniformMDomain<DimVx, DimX>;
using DBlockVxX = Block<MDomainVxX, double>;

class DBlockXTest : public ::testing::Test
{
protected:
    MeshX mesh = MeshX(0.0, 1.0);
    MDomainX dom = MDomainX(mesh, 0, 100);
};

TEST_F(DBlockXTest, constructor)
{
    DBlockX block(dom);
}

TEST_F(DBlockXTest, domain)
{
    DBlockX block(dom);
    ASSERT_EQ(dom, block.domain());
}

TEST_F(DBlockXTest, domainX)
{
    DBlockX block(dom);
    ASSERT_EQ(dom, block.domain<DimX>());
}

TEST_F(DBlockXTest, get_domainX)
{
    DBlockX block(dom);
    ASSERT_EQ(dom, get_domain<DimX>(block));
}

TEST_F(DBlockXTest, access)
{
    DBlockX block(dom);
    for (auto&& ii : block.domain()) {
        ASSERT_EQ(block(ii), block(std::array {ii}));
    }
}

TEST_F(DBlockXTest, deepcopy)
{
    DBlockX block(dom);
    for (auto&& ii : block.domain()) {
        block(ii) = 1.001 * ii;
    }
    DBlockX block2(block.domain());
    deepcopy(block2, block);
    for (auto&& ii : block.domain()) {
        // we expect complete equality, not ASSERT_DOUBLE_EQ: these are copy
        ASSERT_EQ(block2(ii), block(ii));
    }
}

class DBlockXVxTest : public ::testing::Test
{
protected:
    MeshXVx mesh = MeshXVx(RCoordXVx(0.0, 0.0), RCoordXVx(1.0, 1.0));
    MDomainXVx dom = MDomainXVx(mesh, MCoordXVx(0, 0), MCoordXVx(100, 100));
};

TEST_F(DBlockXVxTest, deepcopy)
{
    DBlockXVx block(dom);
    for (auto&& ii : block.domain<DimX>()) {
        for (auto&& jj : block.domain<DimVx>()) {
            block(ii, jj) = 1. * ii + .001 * jj;
        }
    }
    DBlockXVx block2(block.domain());
    deepcopy(block2, block);
    for (auto&& ii : block.domain<DimX>()) {
        for (auto&& jj : block.domain<DimVx>()) {
            // we expect complete equality, not ASSERT_DOUBLE_EQ: these are copy
            ASSERT_EQ(block2(ii, jj), block(ii, jj));
        }
    }
}

TEST_F(DBlockXVxTest, reordering)
{
    DBlockXVx block(dom);
    for (auto&& ii : block.domain<DimX>()) {
        for (auto&& jj : block.domain<DimVx>()) {
            block(ii, jj) = 1. * ii + .001 * jj;
        }
    }

    MeshVxX mesh_reordered = MeshVxX(RCoordVxX(0.0, 0.0), RCoordVxX(1.0, 1.0));
    MDomainVxX dom_reordered = MDomainVxX(mesh_reordered, MCoordVxX(0, 0), MCoordVxX(100, 100));
    DBlockVxX block_reordered(dom_reordered);
    deepcopy(block_reordered, block);
    for (auto&& ii : block.domain<DimX>()) {
        for (auto&& jj : block.domain<DimVx>()) {
            // we expect complete equality, not ASSERT_DOUBLE_EQ: these are copy
            ASSERT_EQ(block_reordered(jj, ii), block(ii, jj));
        }
    }
}

TEST_F(DBlockXVxTest, slice)
{
    DBlockXVx block(dom);
    for (auto&& ii : block.domain<DimX>()) {
        for (auto&& jj : block.domain<DimVx>()) {
            block(ii, jj) = 1. * ii + .001 * jj;
        }
    }
    {
        const DBlockXVx& constref_block = block;
        constexpr auto SLICE_VAL = 1;
        auto&& block_x = constref_block.subblockview(all, SLICE_VAL);
        for (auto&& ii : constref_block.domain<DimX>()) {
            // we expect complete equality, not ASSERT_DOUBLE_EQ: these are copy
            ASSERT_EQ(block_x(ii), constref_block(ii, SLICE_VAL));
        }
        auto&& block_v = constref_block.subblockview(SLICE_VAL, all);
        for (auto&& ii : constref_block.domain<DimVx>()) {
            // we expect complete equality, not ASSERT_DOUBLE_EQ: these are copy
            ASSERT_EQ(block_v(ii), constref_block(SLICE_VAL, ii));
        }
    }
}

TEST_F(DBlockXVxTest, view)
{
    DBlockXVx block(dom);
    for (auto&& ii : block.domain<DimX>()) {
        for (auto&& jj : block.domain<DimVx>()) {
            block(ii, jj) = 1. * ii + .001 * jj;
        }
    }
    auto cview = block.cview();
    for (auto&& ii : block.domain<DimX>()) {
        for (auto&& jj : block.domain<DimVx>()) {
            // we expect complete equality, not ASSERT_DOUBLE_EQ: these are copy
            ASSERT_EQ(cview(ii, jj), block(ii, jj));
        }
    }
}

TEST_F(DBlockXVxTest, tag_rank)
{
    DBlockXVx block(dom);
    ASSERT_EQ(block.tag_rank<DimX>(), 0);
    ASSERT_EQ(block.tag_rank<DimVx>(), 1);
}
