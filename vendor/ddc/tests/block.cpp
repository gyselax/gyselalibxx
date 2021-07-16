#include <iosfwd>
#include <memory>
#include <utility>

#include <gtest/gtest.h>

#include "gtest/gtest_pred_impl.h"

#include "block.h"
#include "blockview.h"
#include "deepcopy.h"
#include "mcoord.h"
#include "mdomain.h"
#include "product_mdomain.h"
#include "product_mesh.h"
#include "rcoord.h"
#include "taggedarray.h"
#include "uniform_mesh.h"

#include <experimental/mdspan>

using namespace std;
using namespace std::experimental;

class DimX;
class DimVx;

using MeshX = UniformMesh<DimX>;
using RCoordX = MeshX::rcoord_type;
using MeshVx = UniformMesh<DimVx>;
using RCoordVx = MeshVx::rcoord_type;
using MDomainX = ProductMDomain<MeshX>;
using DBlockX = Block<MDomainX, double>;
using MeshXVx = ProductMesh<MeshX, MeshVx>;
using RCoordXVx = RCoord<DimX, DimVx>;
using MCoordXVx = MCoord<MeshX, MeshVx>;
using MDomainXVx = ProductMDomain<MeshX, MeshVx>;
using DBlockXVx = Block<MDomainXVx, double>;
using MeshVxX = ProductMesh<MeshVx, MeshX>;
using RCoordVxX = RCoord<DimVx, DimX>;
using MCoordVxX = MCoord<MeshVx, MeshX>;
using MDomainVxX = ProductMDomain<MeshVx, MeshX>;
using DBlockVxX = Block<MDomainVxX, double>;

class DBlockXTest : public ::testing::Test
{
protected:
    MeshX mesh = MeshX(0.0, 1.0);
    MDomainX dom = MDomainX(ProductMesh(mesh), 10, 100);
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
    ASSERT_EQ(dom, block.domain<MeshX>());
}

TEST_F(DBlockXTest, get_domainX)
{
    DBlockX block(dom);
    ASSERT_EQ(dom, get_domain<MeshX>(block));
}

TEST_F(DBlockXTest, access)
{
    DBlockX block(dom);
    for (auto&& ii : block.domain()) {
        ASSERT_EQ(block(ii), block(ii));
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
    MeshX mesh_x = MeshX(0.0, 1.0);
    MeshVx mesh_vx = MeshVx(0.0, 1.0);
    MeshXVx mesh = MeshXVx(mesh_x, mesh_vx);
    MDomainXVx dom = MDomainXVx(mesh, MCoordXVx(0, 0), MCoordXVx(100, 100));
};

TEST_F(DBlockXVxTest, deepcopy)
{
    DBlockXVx block(dom);
    for (auto&& ii : block.domain<MeshX>()) {
        for (auto&& jj : block.domain<MeshVx>()) {
            block(ii, jj) = 1. * ii + .001 * jj;
        }
    }
    DBlockXVx block2(block.domain());
    deepcopy(block2, block);
    for (auto&& ii : block.domain<MeshX>()) {
        for (auto&& jj : block.domain<MeshVx>()) {
            // we expect complete equality, not ASSERT_DOUBLE_EQ: these are copy
            ASSERT_EQ(block2(ii, jj), block(ii, jj));
        }
    }
}

TEST_F(DBlockXVxTest, reordering)
{
    DBlockXVx block(dom);
    for (auto&& ii : block.domain<MeshX>()) {
        for (auto&& jj : block.domain<MeshVx>()) {
            block(ii, jj) = 1. * ii + .001 * jj;
        }
    }

    MDomainVxX dom_reordered = select<MeshVx, MeshX>(dom);
    DBlockVxX block_reordered(dom_reordered);
    deepcopy(block_reordered, block);
    for (auto&& ii : block.domain<MeshX>()) {
        for (auto&& jj : block.domain<MeshVx>()) {
            // we expect complete equality, not ASSERT_DOUBLE_EQ: these are copy
            ASSERT_EQ(block_reordered(jj, ii), block(ii, jj));
        }
    }
}

TEST_F(DBlockXVxTest, slice)
{
    DBlockXVx block(dom);
    for (auto&& ii : block.domain<MeshX>()) {
        for (auto&& jj : block.domain<MeshVx>()) {
            block(ii, jj) = 1. * ii + .001 * jj;
        }
    }
    {
        const DBlockXVx& constref_block = block;
        constexpr auto SLICE_VAL = 1;
        auto&& block_x = constref_block.subblockview(full_extent, SLICE_VAL);
        ASSERT_EQ(block_x.extent(0), block.extent(0));
        for (auto&& ii : constref_block.domain<MeshX>()) {
            // we expect complete equality, not ASSERT_DOUBLE_EQ: these are copy
            ASSERT_EQ(block_x(ii), constref_block(ii, SLICE_VAL));
        }

        auto&& block_v = constref_block.subblockview(SLICE_VAL, full_extent);
        ASSERT_EQ(block_v.extent(0), block.extent(1));
        for (auto&& ii : constref_block.domain<MeshVx>()) {
            // we expect complete equality, not ASSERT_DOUBLE_EQ: these are copy
            ASSERT_EQ(block_v(ii), constref_block(SLICE_VAL, ii));
        }

        auto&& subblock = constref_block.subblockview(std::pair(10, 5), full_extent);
        ASSERT_EQ(subblock.extent(0), 5);
        ASSERT_EQ(subblock.extent(1), get<MeshVx>(block.domain()).size());
        for (auto&& ii : subblock.domain<MeshX>()) {
            for (auto&& jj : subblock.domain<MeshVx>()) {
                // we expect complete equality, not ASSERT_DOUBLE_EQ: these are copy
                ASSERT_EQ(subblock(ii, jj), constref_block(10 + ii, jj));
            }
        }
    }
}

TEST_F(DBlockXVxTest, view)
{
    DBlockXVx block(dom);
    for (auto&& ii : block.domain<MeshX>()) {
        for (auto&& jj : block.domain<MeshVx>()) {
            block(ii, jj) = 1. * ii + .001 * jj;
        }
    }
    auto cview = block.cview();
    for (auto&& ii : block.domain<MeshX>()) {
        for (auto&& jj : block.domain<MeshVx>()) {
            // we expect complete equality, not ASSERT_DOUBLE_EQ: these are copy
            ASSERT_EQ(cview(ii, jj), block(ii, jj));
        }
    }
}
