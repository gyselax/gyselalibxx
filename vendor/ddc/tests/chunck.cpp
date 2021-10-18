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

TEST(Chunk1DTest, layout_type)
{
    ChunkX<double> chunck(ddim_x);

    EXPECT_TRUE((std::is_same_v<
                 std::decay_t<decltype(chunck)>::layout_type,
                 std::experimental::layout_right>));
}

// TODO: many missing types

// \}
// Functions implemented in Chunk 1D (and free functions specific to it) \{

TEST(Chunk1DTest, move_constructor)
{
    double constexpr factor = 1.391;
    ChunkX<double> chunck(ddim_x);
    for (auto&& ix : chunck.domain()) {
        chunck(ix) = factor * ix;
    }

    ChunkX<double> chunck2(std::move(chunck));
    EXPECT_EQ(chunck2.domain(), ddim_x);
    for (auto&& ix : chunck2.domain()) {
        // we expect exact equality, not EXPECT_DOUBLE_EQ: this is the same ref twice
        EXPECT_EQ(factor * ix, chunck2(ix));
    }
}

TEST(Chunk1DTest, move_assignment)
{
    double constexpr factor = 1.976;
    ChunkX<double> chunck(ddim_x);
    for (auto&& ix : chunck.domain()) {
        chunck(ix) = factor * ix;
    }

    ChunkX<double> chunck2(IDomainX(idim_x, lbound_x, IVectX(0)));
    chunck2 = std::move(chunck);
    EXPECT_EQ(chunck2.domain(), ddim_x);
    for (auto&& ix : chunck2.domain()) {
        // we expect exact equality, not EXPECT_DOUBLE_EQ: this is the same ref twice
        EXPECT_EQ(factor * ix, chunck2(ix));
    }
}

TEST(Chunk1DTest, swap)
{
    double constexpr factor = 1.976;
    ChunkX<double> chunck(ddim_x);
    for (auto&& ix : chunck.domain()) {
        chunck(ix) = factor * ix;
    }

    IDomainX empty_domain(idim_x, lbound_x, IVectX(0));
    ChunkX<double> chunck2(empty_domain);

    chunck2.swap(chunck);
    EXPECT_EQ(chunck.domain(), empty_domain);
    EXPECT_EQ(chunck2.domain(), ddim_x);
    for (auto&& ix : chunck2.domain()) {
        // we expect exact equality, not EXPECT_DOUBLE_EQ: this is the same ref twice
        EXPECT_EQ(factor * ix, chunck2(ix));
    }
}

// no dim subset access in 1D

TEST(Chunk1DTest, access_const)
{
    double constexpr factor = 1.012;
    ChunkX<double> chunck(ddim_x);
    ChunkX<double> const& chunck_cref = chunck;
    for (auto&& ix : chunck.domain()) {
        chunck(ix) = factor * ix;
        // we expect exact equality, not EXPECT_DOUBLE_EQ: this is the same ref twice
        EXPECT_EQ(chunck_cref(ix), factor * ix);
    }
}

TEST(Chunk1DTest, access)
{
    double constexpr factor = 1.012;
    ChunkX<double> chunck(ddim_x);
    for (auto&& ix : chunck.domain()) {
        chunck(ix) = factor * ix;
        // we expect exact equality, not EXPECT_DOUBLE_EQ: this is the same ref twice
        EXPECT_EQ(chunck(ix), factor * ix);
    }
}

TEST(Chunk1DTest, span_cview)
{
    double constexpr factor = 1.567;
    ChunkX<double> chunck(ddim_x);
    for (auto&& ix : chunck.domain()) {
        chunck(ix) = factor * ix;
        // we expect exact equality, not EXPECT_DOUBLE_EQ: this is the same ref twice
        EXPECT_EQ(chunck.span_cview()(ix), factor * ix);
    }
}

TEST(Chunk1DTest, view_const)
{
    double constexpr factor = 1.802;
    ChunkX<double> chunck(ddim_x);
    ChunkX<double> const& chunck_cref = chunck;
    for (auto&& ix : chunck.domain()) {
        chunck(ix) = factor * ix;
        // we expect exact equality, not EXPECT_DOUBLE_EQ: this is the same ref twice
        EXPECT_EQ(chunck_cref.span_view()(ix), factor * ix);
    }
}

TEST(Chunk1DTest, view)
{
    double constexpr factor = 1.259;
    ChunkX<double> chunck(ddim_x);
    for (auto&& ix : chunck.domain()) {
        chunck.span_view()(ix) = factor * ix;
        // we expect exact equality, not EXPECT_DOUBLE_EQ: this is the same ref twice
        EXPECT_EQ(chunck(ix), factor * ix);
    }
}

// \}
// Functions inherited from ChunkSpan (and free functions implemented for it) \{

// constructors are hidden

// assignment operators are hidden

// no slicing (operator[]) in 1D

// access (operator()) operators are hidden

TEST(Chunk1DTest, ibegin)
{
    ChunkX<double> chunck(ddim_x);
    EXPECT_EQ(chunck.ibegin<IDimX>(), lbound_x);
}

TEST(Chunk1DTest, iend)
{
    ChunkX<double> chunck(ddim_x);
    EXPECT_EQ(chunck.iend<IDimX>(), sentinel_x);
}

// TODO: accessor

TEST(Chunk1DTest, rank)
{
    ChunkX<double> chunck(ddim_x);
    EXPECT_EQ(ChunkX<double>::rank(), 1);
}

TEST(Chunk1DTest, rank_dynamic)
{
    ChunkX<double> chunck(ddim_x);
    EXPECT_EQ(ChunkX<double>::rank_dynamic(), 1);
}

TEST(Chunk1DTest, extents)
{
    ChunkX<double> chunck(ddim_x);
    EXPECT_EQ(chunck.extents(), npoints_x);
}

TEST(Chunk1DTest, extent)
{
    ChunkX<double> chunck(ddim_x);
    EXPECT_EQ(chunck.extent<IDimX>(), npoints_x.value());
}

TEST(Chunk1DTest, size)
{
    ChunkX<double> chunck(ddim_x);
    EXPECT_EQ(chunck.size(), npoints_x.value());
}

TEST(Chunk1DTest, unique_size)
{
    ChunkX<double> chunck(ddim_x);
    EXPECT_EQ(chunck.unique_size(), npoints_x.value());
}

TEST(Chunk1DTest, is_always_unique)
{
    ChunkX<double> chunck(ddim_x);
    EXPECT_TRUE(chunck.is_always_unique());
}

TEST(Chunk1DTest, is_always_contiguous)
{
    ChunkX<double> chunck(ddim_x);
    EXPECT_TRUE(chunck.is_always_contiguous());
}

TEST(Chunk1DTest, is_always_strided)
{
    ChunkX<double> chunck(ddim_x);
    EXPECT_TRUE(chunck.is_always_strided());
}

// TODO: mapping

TEST(Chunk1DTest, is_unique)
{
    ChunkX<double> chunck(ddim_x);
    EXPECT_TRUE(chunck.is_unique());
}

TEST(Chunk1DTest, is_contiguous)
{
    ChunkX<double> chunck(ddim_x);
    EXPECT_TRUE(chunck.is_contiguous());
}

TEST(Chunk1DTest, is_strided)
{
    ChunkX<double> chunck(ddim_x);
    EXPECT_TRUE(chunck.is_strided());
}

// TODO: stride

// swap is hidden

TEST(Chunk1DTest, domain)
{
    ChunkX<double> chunck(ddim_x);
    EXPECT_EQ(ddim_x, chunck.domain());
}

TEST(Chunk1DTest, domainX)
{
    ChunkX<double> chunck(ddim_x);
    EXPECT_EQ(ddim_x, chunck.domain<IDimX>());
}

// TODO: data

// TODO: internal_mdspan

// TODO: allocation_mdspan

TEST(Chunk1DTest, get_domainX)
{
    ChunkX<double> chunck(ddim_x);
    EXPECT_EQ(ddim_x, get_domain<IDimX>(chunck));
}

TEST(Chunk1DTest, deepcopy)
{
    ChunkX<double> chunck(ddim_x);
    for (auto&& ix : chunck.domain()) {
        chunck(ix) = 1.001 * ix;
    }
    ChunkX<double> chunck2(chunck.domain());
    deepcopy(chunck2, chunck);
    for (auto&& ix : chunck.domain()) {
        // we expect exact equality, not EXPECT_DOUBLE_EQ: these are copy
        EXPECT_EQ(chunck2(ix), chunck(ix));
    }
}

// \}
// Functions implemented in Chunk 2D (and free functions specific to it) \{

// TODO: lots to do still!

TEST(Chunk2DTest, access)
{
    ChunkXY<double> chunck(ddim_x_y);
    for (auto&& ix : chunck.domain<IDimX>()) {
        for (auto&& iy : chunck.domain<IDimY>()) {
            chunck(ix, iy) = 1.357 * ix + 1.159 * iy;
            // we expect exact equality, not EXPECT_DOUBLE_EQ: this is the same ref twice
            EXPECT_EQ(chunck(ix, iy), chunck(ix, iy));
        }
    }
}

TEST(Chunk2DTest, access_reordered)
{
    ChunkXY<double> chunck(ddim_x_y);
    for (auto&& ix : chunck.domain<IDimX>()) {
        for (auto&& iy : chunck.domain<IDimY>()) {
            chunck(ix, iy) = 1.455 * ix + 1.522 * iy;
            // we expect exact equality, not EXPECT_DOUBLE_EQ: this is the same ref twice
            EXPECT_EQ(chunck(iy, ix), chunck(ix, iy));
        }
    }
}

TEST(Chunk2DTest, cview)
{
    ChunkXY<double> chunck(ddim_x_y);
    for (auto&& ix : chunck.domain<IDimX>()) {
        for (auto&& iy : chunck.domain<IDimY>()) {
            chunck(ix, iy) = 1. * ix + .001 * iy;
        }
    }
    auto cview = chunck.span_cview();
    for (auto&& ix : chunck.domain<IDimX>()) {
        for (auto&& iy : chunck.domain<IDimY>()) {
            // we expect complete equality, not EXPECT_DOUBLE_EQ: these are copy
            EXPECT_EQ(cview(ix, iy), chunck(ix, iy));
        }
    }
}

TEST(Chunk2DTest, slice_coord_x)
{
    IndexX constexpr slice_x_val = IndexX(lbound_x + 1);

    ChunkXY<double> chunck(ddim_x_y);
    ChunkXY<double> const& chunck_cref = chunck;
    for (auto&& ix : chunck.domain<IDimX>()) {
        for (auto&& iy : chunck.domain<IDimY>()) {
            chunck(ix, iy) = 1. * ix + .001 * iy;
        }
    }

    auto&& chunck_y = chunck_cref[slice_x_val];
    EXPECT_TRUE((std::is_same_v<
                 std::decay_t<decltype(chunck_y)>::layout_type,
                 std::experimental::layout_right>));
    EXPECT_EQ(chunck_y.extent<IDimY>(), chunck.extent<IDimY>());
    for (auto&& ix : chunck_cref.domain<IDimY>()) {
        // we expect complete equality, not EXPECT_DOUBLE_EQ: these are copy
        EXPECT_EQ(chunck_y(ix), chunck_cref(slice_x_val, ix));
    }
}

TEST(Chunk2DTest, slice_coord_y)
{
    MCoordY constexpr slice_y_val = MCoordY(lbound_y + 1);

    ChunkXY<double> chunck(ddim_x_y);
    ChunkXY<double> const& chunck_cref = chunck;
    for (auto&& ix : chunck.domain<IDimX>()) {
        for (auto&& iy : chunck.domain<IDimY>()) {
            chunck(ix, iy) = 1. * ix + .001 * iy;
        }
    }

    auto&& chunck_x = chunck_cref[slice_y_val];
    EXPECT_TRUE((std::is_same_v<
                 std::decay_t<decltype(chunck_x)>::layout_type,
                 std::experimental::layout_stride>));
    EXPECT_EQ(chunck_x.extent<IDimX>(), chunck.extent<IDimX>());
    for (auto&& ix : chunck_cref.domain<IDimX>()) {
        // we expect complete equality, not EXPECT_DOUBLE_EQ: these are copy
        EXPECT_EQ(chunck_x(ix), chunck_cref(ix, slice_y_val));
    }
}

TEST(Chunk2DTest, slice_domain_x)
{
    IDomainX constexpr subdomain_x = IDomainX(idim_x, IndexX(lbound_x + 1), IVectX(npoints_x - 2));

    ChunkXY<double> chunck(ddim_x_y);
    ChunkXY<double> const& chunck_cref = chunck;
    for (auto&& ix : chunck.domain<IDimX>()) {
        for (auto&& iy : chunck.domain<IDimY>()) {
            chunck(ix, iy) = 1. * ix + .001 * iy;
        }
    }

    auto&& subchunck_x = chunck_cref[subdomain_x];
    EXPECT_TRUE((std::is_same_v<
                 std::decay_t<decltype(subchunck_x)>::layout_type,
                 std::experimental::layout_right>));

    EXPECT_EQ(subchunck_x.extent<IDimX>(), subdomain_x.size());
    EXPECT_EQ(subchunck_x.extent<IDimY>(), chunck.domain<IDimY>().size());
    for (auto&& ix : subchunck_x.domain<IDimX>()) {
        for (auto&& iy : subchunck_x.domain<IDimY>()) {
            // we expect complete equality, not EXPECT_DOUBLE_EQ: these are copy
            EXPECT_EQ(subchunck_x(ix, iy), chunck_cref(ix, iy));
        }
    }
}

TEST(Chunk2DTest, slice_domain_x_tooearly)
{
    IDomainX constexpr subdomain_x = IDomainX(idim_x, IndexX(lbound_x - 1), npoints_x);

    ChunkXY<double> chunck(ddim_x_y);
#ifndef NDEBUG // The assertion is only checked if NDEBUG isn't defined
    // the error message is checked with clang & gcc only
    ASSERT_DEATH(
            chunck[subdomain_x],
            R"rgx([Aa]ssert.*get<ODDims>\(m_lbound\).*get<ODDims>\(odomain\.m_lbound\))rgx");
#endif
}

TEST(Chunk2DTest, slice_domain_x_toolate)
{
    IDomainX constexpr subdomain_x = IDomainX(idim_x, lbound_x, IVectX(npoints_x + 1));

    ChunkXY<double> chunck(ddim_x_y);
#ifndef NDEBUG // The assertion is only checked if NDEBUG isn't defined
    // the error message is checked with clang & gcc only
    ASSERT_DEATH(
            chunck[subdomain_x],
            R"rgx([Aa]ssert.*get<ODDims>\(m_ubound\).*get<ODDims>\(odomain\.m_ubound\).*)rgx");
#endif
}

TEST(Chunk2DTest, slice_domain_y)
{
    MDomainY constexpr subdomain_y
            = MDomainY(idim_y, MCoordY(lbound_y + 1), MLengthY(npoints_y - 2));

    ChunkXY<double> chunck(ddim_x_y);
    ChunkXY<double> const& chunck_cref = chunck;
    for (auto&& ix : chunck.domain<IDimX>()) {
        for (auto&& iy : chunck.domain<IDimY>()) {
            chunck(ix, iy) = 1. * ix + .001 * iy;
        }
    }
    auto&& subchunck_y = chunck_cref[subdomain_y];
    EXPECT_TRUE((std::is_same_v<
                 std::decay_t<decltype(subchunck_y)>::layout_type,
                 std::experimental::layout_stride>));

    EXPECT_EQ(subchunck_y.extent<IDimX>(), chunck.domain<IDimX>().size());
    EXPECT_EQ(subchunck_y.extent<IDimY>(), subdomain_y.size());
    for (auto&& ix : subchunck_y.domain<IDimX>()) {
        for (auto&& iy : subchunck_y.domain<IDimY>()) {
            // we expect complete equality, not EXPECT_DOUBLE_EQ: these are copy
            EXPECT_EQ(subchunck_y(ix, iy), chunck_cref(ix, iy));
        }
    }
}

TEST(Chunk2DTest, deepcopy)
{
    ChunkXY<double> chunck(ddim_x_y);
    for (auto&& ix : chunck.domain<IDimX>()) {
        for (auto&& iy : chunck.domain<IDimY>()) {
            chunck(ix, iy) = 1.739 * ix + 1.412 * iy;
        }
    }
    ChunkXY<double> chunck2(chunck.domain());
    deepcopy(chunck2, chunck);
    for (auto&& ix : chunck.domain<IDimX>()) {
        for (auto&& iy : chunck.domain<IDimY>()) {
            // we expect complete equality, not EXPECT_DOUBLE_EQ: these are copy
            EXPECT_EQ(chunck2(ix, iy), chunck(ix, iy));
        }
    }
}

TEST(Chunk2DTest, deepcopy_reordered)
{
    ChunkXY<double> chunck(ddim_x_y);
    for (auto&& ix : chunck.domain<IDimX>()) {
        for (auto&& iy : chunck.domain<IDimY>()) {
            chunck(ix, iy) = 1.739 * ix + 1.412 * iy;
        }
    }
    ChunkYX<double> chunck2(select<IDimY, IDimX>(chunck.domain()));
    deepcopy(chunck2, chunck);
    for (auto&& ix : chunck.domain<IDimX>()) {
        for (auto&& iy : chunck.domain<IDimY>()) {
            // we expect complete equality, not EXPECT_DOUBLE_EQ: these are copy
            EXPECT_EQ(chunck2(ix, iy), chunck(ix, iy));
            EXPECT_EQ(chunck2(ix, iy), chunck(iy, ix));
        }
    }
}
