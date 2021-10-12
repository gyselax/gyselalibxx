// SPDX-License-Identifier: MIT
#include <iosfwd>
#include <memory>
#include <type_traits>

#include <experimental/mdspan>

#include <ddc/Block>
#include <ddc/BlockSpan>
#include <ddc/MCoord>
#include <ddc/ProductMDomain>
#include <ddc/RCoord>
#include <ddc/TaggedVector>
#include <ddc/UniformMesh>

#include <gtest/gtest.h>

using namespace std;
using namespace std::experimental;

namespace {

class RDimX;
using RCoordX = RCoord<RDimX>;

using IDimX = UniformMesh<RDimX>;
using MCoordX = MCoord<IDimX>;
using MLengthX = MLength<IDimX>;
using MDomainX = ProductMDomain<IDimX>;

template <class Datatype>
using BlockX = Block<Datatype, MDomainX>;
template <class Datatype>
using BlockSpanX = BlockSpan<Datatype, MDomainX>;


class RDimY;
using RCoordY = RCoord<RDimY>;

using IDimY = UniformMesh<RDimY>;
using MCoordY = MCoord<IDimY>;
using MLengthY = MLength<IDimY>;
using MDomainY = ProductMDomain<IDimY>;

template <class Datatype>
using BlockY = Block<Datatype, MDomainY>;


class RDimZ;
using RCoordZ = RCoord<RDimZ>;

using IDimZ = UniformMesh<RDimZ>;
using MCoordZ = MCoord<IDimZ>;
using MLengthZ = MLength<IDimZ>;
using MDomainZ = ProductMDomain<IDimZ>;


using RCoordXY = RCoord<RDimY, RDimX>;

using MCoordXY = MCoord<IDimX, IDimY>;
using MLengthXY = MLength<IDimX, IDimY>;
using MDomainXY = ProductMDomain<IDimX, IDimY>;

template <class Datatype>
using BlockXY = Block<Datatype, MDomainXY>;


using RCoordYX = RCoord<RDimX, RDimY>;

using MCoordYX = MCoord<IDimY, IDimX>;
using MLengthYX = MLength<IDimY, IDimX>;
using MDomainYX = ProductMDomain<IDimY, IDimX>;

template <class Datatype>
using BlockYX = Block<Datatype, MDomainYX>;


static RCoordX constexpr origin_x(0);
static RCoordX constexpr step_x(.01);
static IDimX constexpr idim_x = IDimX(origin_x, step_x);
static MCoordX constexpr lbound_x(50);
static MLengthX constexpr npoints_x(3);
static MCoordX constexpr sentinel_x(lbound_x + npoints_x);
static MDomainX constexpr mesh_x = MDomainX(idim_x, lbound_x, npoints_x);

static RCoordY constexpr origin_y(.1);
static RCoordY constexpr step_y(.103);
static IDimY constexpr idim_y = IDimY(origin_y, step_y);
static MCoordY constexpr lbound_y(4);
static MLengthY constexpr npoints_y(12);

static RCoordZ constexpr origin_z(1.);
static RCoordZ constexpr step_z(.201);
static IDimZ constexpr idim_z = IDimZ(origin_z, step_z);
static MCoordZ constexpr lbound_z(2);
static MLengthZ constexpr npoints_z(10);

static MCoordXY constexpr lbound_x_y {lbound_x, lbound_y};
static MLengthXY constexpr npoints_x_y(npoints_x, npoints_y);
static MDomainXY constexpr mesh_x_y = MDomainXY(idim_x, idim_y, lbound_x_y, npoints_x_y);

} // namespace

// Member types of Block 1D \{

TEST(Block1DTest, layout_type)
{
    BlockX<double> block(mesh_x);

    EXPECT_TRUE((std::is_same_v<
                 std::decay_t<decltype(block)>::layout_type,
                 std::experimental::layout_right>));
}

// TODO: many missing types

// \}
// Functions implemented in Block 1D (and free functions specific to it) \{

TEST(Block1DTest, move_constructor)
{
    double constexpr factor = 1.391;
    BlockX<double> block(mesh_x);
    for (auto&& ix : block.domain()) {
        block(ix) = factor * ix;
    }

    BlockX<double> block2(std::move(block));
    EXPECT_EQ(block2.domain(), mesh_x);
    for (auto&& ix : block2.domain()) {
        // we expect exact equality, not EXPECT_DOUBLE_EQ: this is the same ref twice
        EXPECT_EQ(factor * ix, block2(ix));
    }
}

TEST(Block1DTest, move_assignment)
{
    double constexpr factor = 1.976;
    BlockX<double> block(mesh_x);
    for (auto&& ix : block.domain()) {
        block(ix) = factor * ix;
    }

    BlockX<double> block2(MDomainX(idim_x, lbound_x, MLengthX(0)));
    block2 = std::move(block);
    EXPECT_EQ(block2.domain(), mesh_x);
    for (auto&& ix : block2.domain()) {
        // we expect exact equality, not EXPECT_DOUBLE_EQ: this is the same ref twice
        EXPECT_EQ(factor * ix, block2(ix));
    }
}

TEST(Block1DTest, swap)
{
    double constexpr factor = 1.976;
    BlockX<double> block(mesh_x);
    for (auto&& ix : block.domain()) {
        block(ix) = factor * ix;
    }

    MDomainX empty_domain(idim_x, lbound_x, MLengthX(0));
    BlockX<double> block2(empty_domain);

    block2.swap(block);
    EXPECT_EQ(block.domain(), empty_domain);
    EXPECT_EQ(block2.domain(), mesh_x);
    for (auto&& ix : block2.domain()) {
        // we expect exact equality, not EXPECT_DOUBLE_EQ: this is the same ref twice
        EXPECT_EQ(factor * ix, block2(ix));
    }
}

// no dim subset access in 1D

TEST(Block1DTest, access_const)
{
    double constexpr factor = 1.012;
    BlockX<double> block(mesh_x);
    BlockX<double> const& block_constref = block;
    for (auto&& ix : block.domain()) {
        block(ix) = factor * ix;
        // we expect exact equality, not EXPECT_DOUBLE_EQ: this is the same ref twice
        EXPECT_EQ(block_constref(ix), factor * ix);
    }
}

TEST(Block1DTest, access)
{
    double constexpr factor = 1.012;
    BlockX<double> block(mesh_x);
    for (auto&& ix : block.domain()) {
        block(ix) = factor * ix;
        // we expect exact equality, not EXPECT_DOUBLE_EQ: this is the same ref twice
        EXPECT_EQ(block(ix), factor * ix);
    }
}

TEST(Block1DTest, cview)
{
    double constexpr factor = 1.567;
    BlockX<double> block(mesh_x);
    for (auto&& ix : block.domain()) {
        block(ix) = factor * ix;
        // we expect exact equality, not EXPECT_DOUBLE_EQ: this is the same ref twice
        EXPECT_EQ(block.cview()(ix), factor * ix);
    }
}

TEST(Block1DTest, view_const)
{
    double constexpr factor = 1.802;
    BlockX<double> block(mesh_x);
    BlockX<double> const& block_constref = block;
    for (auto&& ix : block.domain()) {
        block(ix) = factor * ix;
        // we expect exact equality, not EXPECT_DOUBLE_EQ: this is the same ref twice
        EXPECT_EQ(block_constref.view()(ix), factor * ix);
    }
}

TEST(Block1DTest, view)
{
    double constexpr factor = 1.259;
    BlockX<double> block(mesh_x);
    for (auto&& ix : block.domain()) {
        block.view()(ix) = factor * ix;
        // we expect exact equality, not EXPECT_DOUBLE_EQ: this is the same ref twice
        EXPECT_EQ(block(ix), factor * ix);
    }
}

// \}
// Functions inherited from BlockSpan (and free functions implemented for it) \{

// constructors are hidden

// assignment operators are hidden

// no slicing (operator[]) in 1D

// access (operator()) operators are hidden

TEST(Block1DTest, ibegin)
{
    BlockX<double> block(mesh_x);
    EXPECT_EQ(block.ibegin<IDimX>(), lbound_x);
}

TEST(Block1DTest, iend)
{
    BlockX<double> block(mesh_x);
    EXPECT_EQ(block.iend<IDimX>(), sentinel_x);
}

// TODO: accessor

TEST(Block1DTest, rank)
{
    BlockX<double> block(mesh_x);
    EXPECT_EQ(BlockX<double>::rank(), 1);
}

TEST(Block1DTest, rank_dynamic)
{
    BlockX<double> block(mesh_x);
    EXPECT_EQ(BlockX<double>::rank_dynamic(), 1);
}

TEST(Block1DTest, extents)
{
    BlockX<double> block(mesh_x);
    EXPECT_EQ(block.extents(), npoints_x);
}

TEST(Block1DTest, extent)
{
    BlockX<double> block(mesh_x);
    EXPECT_EQ(block.extent<IDimX>(), npoints_x.value());
}

TEST(Block1DTest, size)
{
    BlockX<double> block(mesh_x);
    EXPECT_EQ(block.size(), npoints_x.value());
}

TEST(Block1DTest, unique_size)
{
    BlockX<double> block(mesh_x);
    EXPECT_EQ(block.unique_size(), npoints_x.value());
}

TEST(Block1DTest, is_always_unique)
{
    BlockX<double> block(mesh_x);
    EXPECT_TRUE(block.is_always_unique());
}

TEST(Block1DTest, is_always_contiguous)
{
    BlockX<double> block(mesh_x);
    EXPECT_TRUE(block.is_always_contiguous());
}

TEST(Block1DTest, is_always_strided)
{
    BlockX<double> block(mesh_x);
    EXPECT_TRUE(block.is_always_strided());
}

// TODO: mapping

TEST(Block1DTest, is_unique)
{
    BlockX<double> block(mesh_x);
    EXPECT_TRUE(block.is_unique());
}

TEST(Block1DTest, is_contiguous)
{
    BlockX<double> block(mesh_x);
    EXPECT_TRUE(block.is_contiguous());
}

TEST(Block1DTest, is_strided)
{
    BlockX<double> block(mesh_x);
    EXPECT_TRUE(block.is_strided());
}

// TODO: stride

// swap is hidden

TEST(Block1DTest, domain)
{
    BlockX<double> block(mesh_x);
    EXPECT_EQ(mesh_x, block.domain());
}

TEST(Block1DTest, domainX)
{
    BlockX<double> block(mesh_x);
    EXPECT_EQ(mesh_x, block.domain<IDimX>());
}

// TODO: data

// TODO: internal_mdspan

// TODO: allocation_mdspan

TEST(Block1DTest, get_domainX)
{
    BlockX<double> block(mesh_x);
    EXPECT_EQ(mesh_x, get_domain<IDimX>(block));
}

TEST(Block1DTest, deepcopy)
{
    BlockX<double> block(mesh_x);
    for (auto&& ix : block.domain()) {
        block(ix) = 1.001 * ix;
    }
    BlockX<double> block2(block.domain());
    deepcopy(block2, block);
    for (auto&& ix : block.domain()) {
        // we expect exact equality, not EXPECT_DOUBLE_EQ: these are copy
        EXPECT_EQ(block2(ix), block(ix));
    }
}

// \}
// Functions implemented in Block 2D (and free functions specific to it) \{

// TODO: lots to do still!

TEST(Block2DTest, access)
{
    BlockXY<double> block(mesh_x_y);
    for (auto&& ix : block.domain<IDimX>()) {
        for (auto&& iy : block.domain<IDimY>()) {
            block(ix, iy) = 1.357 * ix + 1.159 * iy;
            // we expect exact equality, not EXPECT_DOUBLE_EQ: this is the same ref twice
            EXPECT_EQ(block(ix, iy), block(ix, iy));
        }
    }
}

TEST(Block2DTest, access_reordered)
{
    BlockXY<double> block(mesh_x_y);
    for (auto&& ix : block.domain<IDimX>()) {
        for (auto&& iy : block.domain<IDimY>()) {
            block(ix, iy) = 1.455 * ix + 1.522 * iy;
            // we expect exact equality, not EXPECT_DOUBLE_EQ: this is the same ref twice
            EXPECT_EQ(block(iy, ix), block(ix, iy));
        }
    }
}

TEST(Block2DTest, cview)
{
    BlockXY<double> block(mesh_x_y);
    for (auto&& ix : block.domain<IDimX>()) {
        for (auto&& iy : block.domain<IDimY>()) {
            block(ix, iy) = 1. * ix + .001 * iy;
        }
    }
    auto cview = block.cview();
    for (auto&& ix : block.domain<IDimX>()) {
        for (auto&& iy : block.domain<IDimY>()) {
            // we expect complete equality, not EXPECT_DOUBLE_EQ: these are copy
            EXPECT_EQ(cview(ix, iy), block(ix, iy));
        }
    }
}

TEST(Block2DTest, slice_coord_x)
{
    MCoordX constexpr slice_x_val = MCoordX(lbound_x + 1);

    BlockXY<double> block(mesh_x_y);
    BlockXY<double> const& block_constref = block;
    for (auto&& ix : block.domain<IDimX>()) {
        for (auto&& iy : block.domain<IDimY>()) {
            block(ix, iy) = 1. * ix + .001 * iy;
        }
    }

    auto&& block_y = block_constref[slice_x_val];
    EXPECT_TRUE((std::is_same_v<
                 std::decay_t<decltype(block_y)>::layout_type,
                 std::experimental::layout_right>));
    EXPECT_EQ(block_y.extent<IDimY>(), block.extent<IDimY>());
    for (auto&& ix : block_constref.domain<IDimY>()) {
        // we expect complete equality, not EXPECT_DOUBLE_EQ: these are copy
        EXPECT_EQ(block_y(ix), block_constref(slice_x_val, ix));
    }
}

TEST(Block2DTest, slice_coord_y)
{
    MCoordY constexpr slice_y_val = MCoordY(lbound_y + 1);

    BlockXY<double> block(mesh_x_y);
    BlockXY<double> const& block_constref = block;
    for (auto&& ix : block.domain<IDimX>()) {
        for (auto&& iy : block.domain<IDimY>()) {
            block(ix, iy) = 1. * ix + .001 * iy;
        }
    }

    auto&& block_x = block_constref[slice_y_val];
    EXPECT_TRUE((std::is_same_v<
                 std::decay_t<decltype(block_x)>::layout_type,
                 std::experimental::layout_stride>));
    EXPECT_EQ(block_x.extent<IDimX>(), block.extent<IDimX>());
    for (auto&& ix : block_constref.domain<IDimX>()) {
        // we expect complete equality, not EXPECT_DOUBLE_EQ: these are copy
        EXPECT_EQ(block_x(ix), block_constref(ix, slice_y_val));
    }
}

TEST(Block2DTest, slice_domain_x)
{
    MDomainX constexpr subdomain_x
            = MDomainX(idim_x, MCoordX(lbound_x + 1), MLengthX(npoints_x - 2));

    BlockXY<double> block(mesh_x_y);
    BlockXY<double> const& block_constref = block;
    for (auto&& ix : block.domain<IDimX>()) {
        for (auto&& iy : block.domain<IDimY>()) {
            block(ix, iy) = 1. * ix + .001 * iy;
        }
    }

    auto&& subblock_x = block_constref[subdomain_x];
    EXPECT_TRUE((std::is_same_v<
                 std::decay_t<decltype(subblock_x)>::layout_type,
                 std::experimental::layout_right>));

    EXPECT_EQ(subblock_x.extent<IDimX>(), subdomain_x.size());
    EXPECT_EQ(subblock_x.extent<IDimY>(), block.domain<IDimY>().size());
    for (auto&& ix : subblock_x.domain<IDimX>()) {
        for (auto&& iy : subblock_x.domain<IDimY>()) {
            // we expect complete equality, not EXPECT_DOUBLE_EQ: these are copy
            EXPECT_EQ(subblock_x(ix, iy), block_constref(ix, iy));
        }
    }
}

TEST(Block2DTest, slice_domain_x_tooearly)
{
    MDomainX constexpr subdomain_x = MDomainX(idim_x, MCoordX(lbound_x - 1), npoints_x);

    BlockXY<double> block(mesh_x_y);
#ifndef NDEBUG // The assertion is only checked if NDEBUG isn't defined
    // the error message is checked with clang & gcc only
    ASSERT_DEATH(
            block[subdomain_x],
            R"rgx([Aa]ssert.*get<OMeshes>\(m_lbound\).*get<OMeshes>\(odomain\.m_lbound\))rgx");
#endif
}

TEST(Block2DTest, slice_domain_x_toolate)
{
    MDomainX constexpr subdomain_x = MDomainX(idim_x, lbound_x, MLengthX(npoints_x + 1));

    BlockXY<double> block(mesh_x_y);
#ifndef NDEBUG // The assertion is only checked if NDEBUG isn't defined
    // the error message is checked with clang & gcc only
    ASSERT_DEATH(
            block[subdomain_x],
            R"rgx([Aa]ssert.*get<OMeshes>\(m_ubound\).*get<OMeshes>\(odomain\.m_ubound\).*)rgx");
#endif
}

TEST(Block2DTest, slice_domain_y)
{
    MDomainY constexpr subdomain_y
            = MDomainY(idim_y, MCoordY(lbound_y + 1), MLengthY(npoints_y - 2));

    BlockXY<double> block(mesh_x_y);
    BlockXY<double> const& block_constref = block;
    for (auto&& ix : block.domain<IDimX>()) {
        for (auto&& iy : block.domain<IDimY>()) {
            block(ix, iy) = 1. * ix + .001 * iy;
        }
    }
    auto&& subblock_y = block_constref[subdomain_y];
    EXPECT_TRUE((std::is_same_v<
                 std::decay_t<decltype(subblock_y)>::layout_type,
                 std::experimental::layout_stride>));

    EXPECT_EQ(subblock_y.extent<IDimX>(), block.domain<IDimX>().size());
    EXPECT_EQ(subblock_y.extent<IDimY>(), subdomain_y.size());
    for (auto&& ix : subblock_y.domain<IDimX>()) {
        for (auto&& iy : subblock_y.domain<IDimY>()) {
            // we expect complete equality, not EXPECT_DOUBLE_EQ: these are copy
            EXPECT_EQ(subblock_y(ix, iy), block_constref(ix, iy));
        }
    }
}

TEST(Block2DTest, deepcopy)
{
    BlockXY<double> block(mesh_x_y);
    for (auto&& ix : block.domain<IDimX>()) {
        for (auto&& iy : block.domain<IDimY>()) {
            block(ix, iy) = 1.739 * ix + 1.412 * iy;
        }
    }
    BlockXY<double> block2(block.domain());
    deepcopy(block2, block);
    for (auto&& ix : block.domain<IDimX>()) {
        for (auto&& iy : block.domain<IDimY>()) {
            // we expect complete equality, not EXPECT_DOUBLE_EQ: these are copy
            EXPECT_EQ(block2(ix, iy), block(ix, iy));
        }
    }
}

TEST(Block2DTest, deepcopy_reordered)
{
    BlockXY<double> block(mesh_x_y);
    for (auto&& ix : block.domain<IDimX>()) {
        for (auto&& iy : block.domain<IDimY>()) {
            block(ix, iy) = 1.739 * ix + 1.412 * iy;
        }
    }
    BlockYX<double> block2(select<IDimY, IDimX>(block.domain()));
    deepcopy(block2, block);
    for (auto&& ix : block.domain<IDimX>()) {
        for (auto&& iy : block.domain<IDimY>()) {
            // we expect complete equality, not EXPECT_DOUBLE_EQ: these are copy
            EXPECT_EQ(block2(ix, iy), block(ix, iy));
            EXPECT_EQ(block2(ix, iy), block(iy, ix));
        }
    }
}
