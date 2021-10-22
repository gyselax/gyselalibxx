// SPDX-License-Identifier: MIT
#include <iosfwd>
#include <memory>
#include <type_traits>

#include <experimental/mdspan>

#include <ddc/Chunk>
#include <ddc/ChunkSpan>
#include <ddc/Coordinate>
#include <ddc/DiscreteCoordinate>
#include <ddc/DiscreteDomain>
#include <ddc/UniformDiscretization>

#include <gtest/gtest.h>

using namespace std;
using namespace std::experimental;

namespace {

class RDimX;
using CoordX = Coordinate<RDimX>;

using IDimX = UniformDiscretization<RDimX>;
using IndexX = DiscreteCoordinate<IDimX>;
using IVectX = DiscreteVector<IDimX>;
using IDomainX = DiscreteDomain<IDimX>;

template <class Datatype>
using ChunkX = Chunk<Datatype, IDomainX>;
template <class Datatype>
using ChunkSpanX = ChunkSpan<Datatype, IDomainX>;


class RDimY;
using RCoordY = Coordinate<RDimY>;

using IDimY = UniformDiscretization<RDimY>;
using MCoordY = DiscreteCoordinate<IDimY>;
using MLengthY = DiscreteVector<IDimY>;
using MDomainY = DiscreteDomain<IDimY>;

template <class Datatype>
using ChunkY = Chunk<Datatype, MDomainY>;


class RDimZ;
using RCoordZ = Coordinate<RDimZ>;

using IDimZ = UniformDiscretization<RDimZ>;
using MCoordZ = DiscreteCoordinate<IDimZ>;
using MLengthZ = DiscreteVector<IDimZ>;
using MDomainZ = DiscreteDomain<IDimZ>;


using RCoordXY = Coordinate<RDimY, RDimX>;

using MCoordXY = DiscreteCoordinate<IDimX, IDimY>;
using MLengthXY = DiscreteVector<IDimX, IDimY>;
using MDomainXY = DiscreteDomain<IDimX, IDimY>;

template <class Datatype>
using ChunkXY = Chunk<Datatype, MDomainXY>;


using RCoordYX = Coordinate<RDimX, RDimY>;

using MCoordYX = DiscreteCoordinate<IDimY, IDimX>;
using MLengthYX = DiscreteVector<IDimY, IDimX>;
using MDomainYX = DiscreteDomain<IDimY, IDimX>;

template <class Datatype>
using ChunkYX = Chunk<Datatype, MDomainYX>;


static CoordX constexpr origin_x(0);
static CoordX constexpr step_x(.01);
static IDimX constexpr idim_x = IDimX(origin_x, step_x);
static IndexX constexpr lbound_x(50);
static IVectX constexpr npoints_x(3);
static IndexX constexpr sentinel_x(lbound_x + npoints_x);
static IDomainX constexpr ddim_x = IDomainX(idim_x, lbound_x, npoints_x);

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
static MDomainXY constexpr ddim_x_y = MDomainXY(idim_x, idim_y, lbound_x_y, npoints_x_y);

} // namespace

// Member types of Chunk 1D \{

TEST(Chunk1DTest, LayoutType)
{
    ChunkX<double> chunk(ddim_x);

    EXPECT_TRUE((std::is_same_v<
                 std::decay_t<decltype(chunk)>::layout_type,
                 std::experimental::layout_right>));
}

// TODO: many missing types

// \}
// Functions implemented in Chunk 1D (and free functions specific to it) \{

TEST(Chunk1DTest, MoveConstructor)
{
    double constexpr factor = 1.391;
    ChunkX<double> chunk(ddim_x);
    for (auto&& ix : chunk.domain()) {
        chunk(ix) = factor * ix;
    }

    ChunkX<double> chunk2(std::move(chunk));
    EXPECT_EQ(chunk2.domain(), ddim_x);
    for (auto&& ix : chunk2.domain()) {
        // we expect exact equality, not EXPECT_DOUBLE_EQ: this is the same ref twice
        EXPECT_EQ(factor * ix, chunk2(ix));
    }
}

TEST(Chunk1DTest, MoveAssignment)
{
    double constexpr factor = 1.976;
    ChunkX<double> chunk(ddim_x);
    for (auto&& ix : chunk.domain()) {
        chunk(ix) = factor * ix;
    }

    ChunkX<double> chunk2(IDomainX(idim_x, lbound_x, IVectX(0)));
    chunk2 = std::move(chunk);
    EXPECT_EQ(chunk2.domain(), ddim_x);
    for (auto&& ix : chunk2.domain()) {
        // we expect exact equality, not EXPECT_DOUBLE_EQ: this is the same ref twice
        EXPECT_EQ(factor * ix, chunk2(ix));
    }
}

TEST(Chunk1DTest, Swap)
{
    double constexpr factor = 1.976;
    ChunkX<double> chunk(ddim_x);
    for (auto&& ix : chunk.domain()) {
        chunk(ix) = factor * ix;
    }

    IDomainX empty_domain(idim_x, lbound_x, IVectX(0));
    ChunkX<double> chunk2(empty_domain);

    chunk2.swap(chunk);
    EXPECT_EQ(chunk.domain(), empty_domain);
    EXPECT_EQ(chunk2.domain(), ddim_x);
    for (auto&& ix : chunk2.domain()) {
        // we expect exact equality, not EXPECT_DOUBLE_EQ: this is the same ref twice
        EXPECT_EQ(factor * ix, chunk2(ix));
    }
}

// no dim subset access in 1D

TEST(Chunk1DTest, AccessConst)
{
    double constexpr factor = 1.012;
    ChunkX<double> chunk(ddim_x);
    ChunkX<double> const& chunk_cref = chunk;
    for (auto&& ix : chunk.domain()) {
        chunk(ix) = factor * ix;
        // we expect exact equality, not EXPECT_DOUBLE_EQ: this is the same ref twice
        EXPECT_EQ(chunk_cref(ix), factor * ix);
    }
}

TEST(Chunk1DTest, Access)
{
    double constexpr factor = 1.012;
    ChunkX<double> chunk(ddim_x);
    for (auto&& ix : chunk.domain()) {
        chunk(ix) = factor * ix;
        // we expect exact equality, not EXPECT_DOUBLE_EQ: this is the same ref twice
        EXPECT_EQ(chunk(ix), factor * ix);
    }
}

TEST(Chunk1DTest, SpanCview)
{
    double constexpr factor = 1.567;
    ChunkX<double> chunk(ddim_x);
    for (auto&& ix : chunk.domain()) {
        chunk(ix) = factor * ix;
        // we expect exact equality, not EXPECT_DOUBLE_EQ: this is the same ref twice
        EXPECT_EQ(chunk.span_cview()(ix), factor * ix);
    }
}

TEST(Chunk1DTest, ViewConst)
{
    double constexpr factor = 1.802;
    ChunkX<double> chunk(ddim_x);
    ChunkX<double> const& chunk_cref = chunk;
    for (auto&& ix : chunk.domain()) {
        chunk(ix) = factor * ix;
        // we expect exact equality, not EXPECT_DOUBLE_EQ: this is the same ref twice
        EXPECT_EQ(chunk_cref.span_view()(ix), factor * ix);
    }
}

TEST(Chunk1DTest, View)
{
    double constexpr factor = 1.259;
    ChunkX<double> chunk(ddim_x);
    for (auto&& ix : chunk.domain()) {
        chunk.span_view()(ix) = factor * ix;
        // we expect exact equality, not EXPECT_DOUBLE_EQ: this is the same ref twice
        EXPECT_EQ(chunk(ix), factor * ix);
    }
}

// \}
// Functions inherited from ChunkSpan (and free functions implemented for it) \{

// constructors are hidden

// assignment operators are hidden

// no slicing (operator[]) in 1D

// access (operator()) operators are hidden

TEST(Chunk1DTest, Ibegin)
{
    ChunkX<double> chunk(ddim_x);
    EXPECT_EQ(chunk.ibegin<IDimX>(), lbound_x);
}

TEST(Chunk1DTest, Iend)
{
    ChunkX<double> chunk(ddim_x);
    EXPECT_EQ(chunk.iend<IDimX>(), sentinel_x);
}

// TODO: accessor

TEST(Chunk1DTest, Rank)
{
    ChunkX<double> chunk(ddim_x);
    EXPECT_EQ(ChunkX<double>::rank(), 1);
}

TEST(Chunk1DTest, RankDynamic)
{
    ChunkX<double> chunk(ddim_x);
    EXPECT_EQ(ChunkX<double>::rank_dynamic(), 1);
}

TEST(Chunk1DTest, Extents)
{
    ChunkX<double> chunk(ddim_x);
    EXPECT_EQ(chunk.extents(), npoints_x);
}

TEST(Chunk1DTest, Extent)
{
    ChunkX<double> chunk(ddim_x);
    EXPECT_EQ(chunk.extent<IDimX>(), npoints_x.value());
}

TEST(Chunk1DTest, Size)
{
    ChunkX<double> chunk(ddim_x);
    EXPECT_EQ(chunk.size(), npoints_x.value());
}

TEST(Chunk1DTest, UniqueSize)
{
    ChunkX<double> chunk(ddim_x);
    EXPECT_EQ(chunk.unique_size(), npoints_x.value());
}

TEST(Chunk1DTest, IsAlwaysUnique)
{
    ChunkX<double> chunk(ddim_x);
    EXPECT_TRUE(chunk.is_always_unique());
}

TEST(Chunk1DTest, IsAlwaysContiguous)
{
    ChunkX<double> chunk(ddim_x);
    EXPECT_TRUE(chunk.is_always_contiguous());
}

TEST(Chunk1DTest, is_always_strided)
{
    ChunkX<double> chunk(ddim_x);
    EXPECT_TRUE(chunk.is_always_strided());
}

// TODO: mapping

TEST(Chunk1DTest, IsUnique)
{
    ChunkX<double> chunk(ddim_x);
    EXPECT_TRUE(chunk.is_unique());
}

TEST(Chunk1DTest, IsContiguous)
{
    ChunkX<double> chunk(ddim_x);
    EXPECT_TRUE(chunk.is_contiguous());
}

TEST(Chunk1DTest, IsStrided)
{
    ChunkX<double> chunk(ddim_x);
    EXPECT_TRUE(chunk.is_strided());
}

// TODO: stride

// swap is hidden

TEST(Chunk1DTest, Domain)
{
    ChunkX<double> chunk(ddim_x);
    EXPECT_EQ(ddim_x, chunk.domain());
}

TEST(Chunk1DTest, DomainX)
{
    ChunkX<double> chunk(ddim_x);
    EXPECT_EQ(ddim_x, chunk.domain<IDimX>());
}

// TODO: data

// TODO: internal_mdspan

// TODO: allocation_mdspan

TEST(Chunk1DTest, GetDomainX)
{
    ChunkX<double> chunk(ddim_x);
    EXPECT_EQ(ddim_x, get_domain<IDimX>(chunk));
}

TEST(Chunk1DTest, Deepcopy)
{
    ChunkX<double> chunk(ddim_x);
    for (auto&& ix : chunk.domain()) {
        chunk(ix) = 1.001 * ix;
    }
    ChunkX<double> chunk2(chunk.domain());
    deepcopy(chunk2, chunk);
    for (auto&& ix : chunk.domain()) {
        // we expect exact equality, not EXPECT_DOUBLE_EQ: these are copy
        EXPECT_EQ(chunk2(ix), chunk(ix));
    }
}

// \}
// Functions implemented in Chunk 2D (and free functions specific to it) \{

// TODO: lots to do still!

TEST(Chunk2DTest, Access)
{
    ChunkXY<double> chunk(ddim_x_y);
    for (auto&& ix : chunk.domain<IDimX>()) {
        for (auto&& iy : chunk.domain<IDimY>()) {
            chunk(ix, iy) = 1.357 * ix + 1.159 * iy;
            // we expect exact equality, not EXPECT_DOUBLE_EQ: this is the same ref twice
            EXPECT_EQ(chunk(ix, iy), chunk(ix, iy));
        }
    }
}

TEST(Chunk2DTest, AccessReordered)
{
    ChunkXY<double> chunk(ddim_x_y);
    for (auto&& ix : chunk.domain<IDimX>()) {
        for (auto&& iy : chunk.domain<IDimY>()) {
            chunk(ix, iy) = 1.455 * ix + 1.522 * iy;
            // we expect exact equality, not EXPECT_DOUBLE_EQ: this is the same ref twice
            EXPECT_EQ(chunk(iy, ix), chunk(ix, iy));
        }
    }
}

TEST(Chunk2DTest, Cview)
{
    ChunkXY<double> chunk(ddim_x_y);
    for (auto&& ix : chunk.domain<IDimX>()) {
        for (auto&& iy : chunk.domain<IDimY>()) {
            chunk(ix, iy) = 1. * ix + .001 * iy;
        }
    }
    auto cview = chunk.span_cview();
    for (auto&& ix : chunk.domain<IDimX>()) {
        for (auto&& iy : chunk.domain<IDimY>()) {
            // we expect complete equality, not EXPECT_DOUBLE_EQ: these are copy
            EXPECT_EQ(cview(ix, iy), chunk(ix, iy));
        }
    }
}

TEST(Chunk2DTest, SliceCoordX)
{
    IndexX constexpr slice_x_val = IndexX(lbound_x + 1);

    ChunkXY<double> chunk(ddim_x_y);
    ChunkXY<double> const& chunk_cref = chunk;
    for (auto&& ix : chunk.domain<IDimX>()) {
        for (auto&& iy : chunk.domain<IDimY>()) {
            chunk(ix, iy) = 1. * ix + .001 * iy;
        }
    }

    auto&& chunk_y = chunk_cref[slice_x_val];
    EXPECT_TRUE((std::is_same_v<
                 std::decay_t<decltype(chunk_y)>::layout_type,
                 std::experimental::layout_right>));
    EXPECT_EQ(chunk_y.extent<IDimY>(), chunk.extent<IDimY>());
    for (auto&& ix : chunk_cref.domain<IDimY>()) {
        // we expect complete equality, not EXPECT_DOUBLE_EQ: these are copy
        EXPECT_EQ(chunk_y(ix), chunk_cref(slice_x_val, ix));
    }
}

TEST(Chunk2DTest, SliceCoordY)
{
    MCoordY constexpr slice_y_val = MCoordY(lbound_y + 1);

    ChunkXY<double> chunk(ddim_x_y);
    ChunkXY<double> const& chunk_cref = chunk;
    for (auto&& ix : chunk.domain<IDimX>()) {
        for (auto&& iy : chunk.domain<IDimY>()) {
            chunk(ix, iy) = 1. * ix + .001 * iy;
        }
    }

    auto&& chunk_x = chunk_cref[slice_y_val];
    EXPECT_TRUE((std::is_same_v<
                 std::decay_t<decltype(chunk_x)>::layout_type,
                 std::experimental::layout_stride>));
    EXPECT_EQ(chunk_x.extent<IDimX>(), chunk.extent<IDimX>());
    for (auto&& ix : chunk_cref.domain<IDimX>()) {
        // we expect complete equality, not EXPECT_DOUBLE_EQ: these are copy
        EXPECT_EQ(chunk_x(ix), chunk_cref(ix, slice_y_val));
    }
}

TEST(Chunk2DTest, SliceDomainX)
{
    IDomainX constexpr subdomain_x = IDomainX(idim_x, IndexX(lbound_x + 1), IVectX(npoints_x - 2));

    ChunkXY<double> chunk(ddim_x_y);
    ChunkXY<double> const& chunk_cref = chunk;
    for (auto&& ix : chunk.domain<IDimX>()) {
        for (auto&& iy : chunk.domain<IDimY>()) {
            chunk(ix, iy) = 1. * ix + .001 * iy;
        }
    }

    auto&& subchunk_x = chunk_cref[subdomain_x];
    EXPECT_TRUE((std::is_same_v<
                 std::decay_t<decltype(subchunk_x)>::layout_type,
                 std::experimental::layout_right>));

    EXPECT_EQ(subchunk_x.extent<IDimX>(), subdomain_x.size());
    EXPECT_EQ(subchunk_x.extent<IDimY>(), chunk.domain<IDimY>().size());
    for (auto&& ix : subchunk_x.domain<IDimX>()) {
        for (auto&& iy : subchunk_x.domain<IDimY>()) {
            // we expect complete equality, not EXPECT_DOUBLE_EQ: these are copy
            EXPECT_EQ(subchunk_x(ix, iy), chunk_cref(ix, iy));
        }
    }
}

TEST(Chunk2DTest, SliceDomainXTooearly)
{
    IDomainX constexpr subdomain_x = IDomainX(idim_x, IndexX(lbound_x - 1), npoints_x);

    ChunkXY<double> chunk(ddim_x_y);
#ifndef NDEBUG // The assertion is only checked if NDEBUG isn't defined
    // the error message is checked with clang & gcc only
    ASSERT_DEATH(
            chunk[subdomain_x],
            R"rgx([Aa]ssert.*get<ODDims>\(m_lbound\).*get<ODDims>\(odomain\.m_lbound\))rgx");
#endif
}

TEST(Chunk2DTest, SliceDomainXToolate)
{
    IDomainX constexpr subdomain_x = IDomainX(idim_x, lbound_x, IVectX(npoints_x + 1));

    ChunkXY<double> chunk(ddim_x_y);
#ifndef NDEBUG // The assertion is only checked if NDEBUG isn't defined
    // the error message is checked with clang & gcc only
    ASSERT_DEATH(
            chunk[subdomain_x],
            R"rgx([Aa]ssert.*get<ODDims>\(m_ubound\).*get<ODDims>\(odomain\.m_ubound\).*)rgx");
#endif
}

TEST(Chunk2DTest, SliceDomainY)
{
    MDomainY constexpr subdomain_y
            = MDomainY(idim_y, MCoordY(lbound_y + 1), MLengthY(npoints_y - 2));

    ChunkXY<double> chunk(ddim_x_y);
    ChunkXY<double> const& chunk_cref = chunk;
    for (auto&& ix : chunk.domain<IDimX>()) {
        for (auto&& iy : chunk.domain<IDimY>()) {
            chunk(ix, iy) = 1. * ix + .001 * iy;
        }
    }
    auto&& subchunk_y = chunk_cref[subdomain_y];
    EXPECT_TRUE((std::is_same_v<
                 std::decay_t<decltype(subchunk_y)>::layout_type,
                 std::experimental::layout_stride>));

    EXPECT_EQ(subchunk_y.extent<IDimX>(), chunk.domain<IDimX>().size());
    EXPECT_EQ(subchunk_y.extent<IDimY>(), subdomain_y.size());
    for (auto&& ix : subchunk_y.domain<IDimX>()) {
        for (auto&& iy : subchunk_y.domain<IDimY>()) {
            // we expect complete equality, not EXPECT_DOUBLE_EQ: these are copy
            EXPECT_EQ(subchunk_y(ix, iy), chunk_cref(ix, iy));
        }
    }
}

TEST(Chunk2DTest, Deepcopy)
{
    ChunkXY<double> chunk(ddim_x_y);
    for (auto&& ix : chunk.domain<IDimX>()) {
        for (auto&& iy : chunk.domain<IDimY>()) {
            chunk(ix, iy) = 1.739 * ix + 1.412 * iy;
        }
    }
    ChunkXY<double> chunk2(chunk.domain());
    deepcopy(chunk2, chunk);
    for (auto&& ix : chunk.domain<IDimX>()) {
        for (auto&& iy : chunk.domain<IDimY>()) {
            // we expect complete equality, not EXPECT_DOUBLE_EQ: these are copy
            EXPECT_EQ(chunk2(ix, iy), chunk(ix, iy));
        }
    }
}

TEST(Chunk2DTest, DeepcopyReordered)
{
    ChunkXY<double> chunk(ddim_x_y);
    for (auto&& ix : chunk.domain<IDimX>()) {
        for (auto&& iy : chunk.domain<IDimY>()) {
            chunk(ix, iy) = 1.739 * ix + 1.412 * iy;
        }
    }
    ChunkYX<double> chunk2(select<IDimY, IDimX>(chunk.domain()));
    deepcopy(chunk2, chunk);
    for (auto&& ix : chunk.domain<IDimX>()) {
        for (auto&& iy : chunk.domain<IDimY>()) {
            // we expect complete equality, not EXPECT_DOUBLE_EQ: these are copy
            EXPECT_EQ(chunk2(ix, iy), chunk(ix, iy));
            EXPECT_EQ(chunk2(ix, iy), chunk(iy, ix));
        }
    }
}
