// SPDX-License-Identifier: MIT
#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include <directional_tag.hpp>
#include <vector_field.hpp>
#include <vector_field_span.hpp>

namespace {

class Tag1
{
};
class Tag2
{
};

using Direction = NDTag<Tag1, Tag2>;

using Coord = ddc::Coordinate<Tag1, Tag2>;

using DElem0D = ddc::DiscreteElement<>;
using DVect0D = ddc::DiscreteVector<>;
using DDom0D = ddc::DiscreteDomain<>;

template <class Datatype>
using VectorField0D = VectorField<Datatype, DDom0D, Direction>;
template <class Datatype>
using VectorFieldSpan0D = VectorFieldSpan<Datatype, DDom0D, Direction>;

struct DDimX;
using DElemX = ddc::DiscreteElement<DDimX>;
using DVectX = ddc::DiscreteVector<DDimX>;
using DDomX = ddc::DiscreteDomain<DDimX>;

template <class Datatype>
using VectorFieldX = VectorField<Datatype, DDomX, Direction>;
template <class Datatype>
using VectorFieldSpanX = VectorFieldSpan<Datatype, DDomX, Direction>;


struct DDimY;
using DElemY = ddc::DiscreteElement<DDimY>;
using DVectY = ddc::DiscreteVector<DDimY>;
using DDomY = ddc::DiscreteDomain<DDimY>;

template <class Datatype>
using VectorFieldY = VectorField<Datatype, DDomY, Direction>;
template <class Datatype>
using VectorFieldSpanY = VectorFieldSpan<Datatype, DDomY, Direction>;


struct DDimZ;
using DElemZ = ddc::DiscreteElement<DDimZ>;
using DVectZ = ddc::DiscreteVector<DDimZ>;
using DDomZ = ddc::DiscreteDomain<DDimZ>;


using DElemXY = ddc::DiscreteElement<DDimX, DDimY>;
using DVectXY = ddc::DiscreteVector<DDimX, DDimY>;
using DDomXY = ddc::DiscreteDomain<DDimX, DDimY>;

template <class Datatype>
using VectorFieldXY = VectorField<Datatype, DDomXY, Direction>;


using DElemYX = ddc::DiscreteElement<DDimY, DDimX>;
using DVectYX = ddc::DiscreteVector<DDimY, DDimX>;
using DDomYX = ddc::DiscreteDomain<DDimY, DDimX>;

template <class Datatype>
using VectorFieldYX = VectorField<Datatype, DDomYX, Direction>;
template <class Datatype>
using VectorFieldSpanYX = VectorFieldSpan<Datatype, DDomYX, Direction>;


static DElem0D constexpr lbound_0d {};
static DVect0D constexpr nelems_0d {};
static DDom0D constexpr dom_0d(lbound_0d, nelems_0d);

static DElemX constexpr lbound_x(50);
static DVectX constexpr nelems_x(3);
static DDomX constexpr dom_x(lbound_x, nelems_x);

static DElemY constexpr lbound_y(4);
static DVectY constexpr nelems_y(12);

static DElemXY constexpr lbound_x_y {lbound_x, lbound_y};
static DVectXY constexpr nelems_x_y(nelems_x, nelems_y);
static DDomXY constexpr dom_x_y(lbound_x_y, nelems_x_y);

} // namespace

// Member types of VectorField 1D \{

TEST(VectorField0DTest, LayoutType)
{
    EXPECT_TRUE((std::is_same_v<
                 VectorField0D<double>::chunk_type::layout_type,
                 std::experimental::layout_right>));
}

TEST(VectorField1DTest, LayoutType)
{
    VectorFieldX<double> chunk(dom_x);

    EXPECT_TRUE((std::is_same_v<
                 std::decay_t<decltype(chunk)::chunk_type>::layout_type,
                 std::experimental::layout_right>));
}

// TODO: many missing types

// \}
// Functions implemented in VectorField 1D (and free functions specific to it) \{

TEST(VectorField0DTest, MoveConstructor)
{
    Coord constexpr factor(1.391, 2.444);
    VectorField0D<double> chunk(dom_0d);
    ddcHelper::get<Tag1>(chunk)() = ddc::get<Tag1>(factor);
    ddcHelper::get<Tag2>(chunk)() = ddc::get<Tag2>(factor);

    VectorField0D<double> chunk2(std::move(chunk));
    EXPECT_EQ(chunk2.domain(), dom_0d);
    EXPECT_DOUBLE_EQ(ddc::get<Tag1>(factor), ddc::get<Tag1>(chunk2()));
    EXPECT_DOUBLE_EQ(ddc::get<Tag2>(factor), ddc::get<Tag2>(chunk2()));
}

TEST(VectorField1DTest, MoveConstructor)
{
    Coord constexpr factor(1.391, 2.444);
    VectorFieldX<double> chunk(dom_x);
    for (auto&& ix : chunk.domain()) {
        Coord val = double(ix.uid()) * factor;
        ddcHelper::get<Tag1>(chunk)(ix) = ddc::get<Tag1>(val);
        ddcHelper::get<Tag2>(chunk)(ix) = ddc::get<Tag2>(val);
    }

    VectorFieldX<double> chunk2(std::move(chunk));
    EXPECT_EQ(chunk2.domain(), dom_x);
    for (auto&& ix : chunk2.domain()) {
        // we expect exact equality, not EXPECT_DOUBLE_EQ: this is the same ref twice
        Coord val = double(ix.uid()) * factor;
        EXPECT_EQ(ddc::get<Tag1>(val), ddc::get<Tag1>(chunk2(ix)));
        EXPECT_EQ(ddc::get<Tag2>(val), ddc::get<Tag2>(chunk2(ix)));
    }
}

TEST(VectorField0DTest, MoveAssignment)
{
    Coord constexpr factor(1.391, 2.444);
    VectorField0D<double> chunk(dom_0d);
    ddcHelper::get<Tag1>(chunk)() = ddc::get<Tag1>(factor);
    ddcHelper::get<Tag2>(chunk)() = ddc::get<Tag2>(factor);

    VectorField0D<double> chunk2(DDom0D(lbound_0d, DVect0D()));
    chunk2 = std::move(chunk);
    EXPECT_EQ(chunk2.domain(), dom_0d);
    EXPECT_DOUBLE_EQ(ddc::get<Tag1>(factor), ddc::get<Tag1>(chunk2()));
    EXPECT_DOUBLE_EQ(ddc::get<Tag2>(factor), ddc::get<Tag2>(chunk2()));
}

TEST(VectorField1DTest, MoveAssignment)
{
    Coord constexpr factor(1.391, 2.444);
    VectorFieldX<double> chunk(dom_x);
    for (auto&& ix : chunk.domain()) {
        Coord val = double(ix.uid()) * factor;
        ddcHelper::get<Tag1>(chunk)(ix) = ddc::get<Tag1>(val);
        ddcHelper::get<Tag2>(chunk)(ix) = ddc::get<Tag2>(val);
    }

    VectorFieldX<double> chunk2(DDomX(lbound_x, DVectX(0)));
    chunk2 = std::move(chunk);
    EXPECT_EQ(chunk2.domain(), dom_x);
    for (auto&& ix : chunk2.domain()) {
        // we expect exact equality, not EXPECT_DOUBLE_EQ: this is the same ref twice
        Coord val = double(ix.uid()) * factor;
        EXPECT_EQ(ddc::get<Tag1>(val), ddc::get<Tag1>(chunk2(ix)));
        EXPECT_EQ(ddc::get<Tag2>(val), ddc::get<Tag2>(chunk2(ix)));
    }
}

TEST(VectorField0DTest, Swap)
{
    Coord constexpr factor(1.391, 2.444);
    VectorField0D<double> chunk(dom_0d);
    ddcHelper::get<Tag1>(chunk)() = ddc::get<Tag1>(factor);
    ddcHelper::get<Tag2>(chunk)() = ddc::get<Tag2>(factor);

    DDom0D empty_domain(lbound_0d, DVect0D());
    VectorField0D<double> chunk2(empty_domain);

    std::swap(chunk2, chunk);
    EXPECT_EQ(chunk.domain(), empty_domain);
    EXPECT_EQ(chunk2.domain(), dom_0d);
    EXPECT_DOUBLE_EQ(ddc::get<Tag1>(factor), ddc::get<Tag1>(chunk2()));
    EXPECT_DOUBLE_EQ(ddc::get<Tag2>(factor), ddc::get<Tag2>(chunk2()));
}

TEST(VectorField1DTest, Swap)
{
    Coord constexpr factor(1.391, 2.444);
    VectorFieldX<double> chunk(dom_x);
    for (auto&& ix : chunk.domain()) {
        Coord val = double(ix.uid()) * factor;
        ddcHelper::get<Tag1>(chunk)(ix) = ddc::get<Tag1>(val);
        ddcHelper::get<Tag2>(chunk)(ix) = ddc::get<Tag2>(val);
    }

    DDomX empty_domain(lbound_x, DVectX(0));
    VectorFieldX<double> chunk2(empty_domain);

    std::swap(chunk2, chunk);
    EXPECT_EQ(chunk.domain(), empty_domain);
    EXPECT_EQ(chunk2.domain(), dom_x);
    for (auto&& ix : chunk2.domain()) {
        // we expect exact equality, not EXPECT_DOUBLE_EQ: this is the same ref twice
        Coord val = double(ix.uid()) * factor;
        EXPECT_EQ(ddc::get<Tag1>(val), ddc::get<Tag1>(chunk2(ix)));
        EXPECT_EQ(ddc::get<Tag2>(val), ddc::get<Tag2>(chunk2(ix)));
    }
}

// no dim subset access in 1D

TEST(VectorField1DTest, AccessConst)
{
    Coord constexpr factor(1.391, 2.444);
    VectorFieldX<double> chunk(dom_x);
    VectorFieldX<double> const& chunk_cref = chunk;
    for (auto&& ix : chunk.domain()) {
        Coord val = double(ix.uid()) * factor;
        ddcHelper::get<Tag1>(chunk)(ix) = ddc::get<Tag1>(val);
        ddcHelper::get<Tag2>(chunk)(ix) = ddc::get<Tag2>(val);
        // we expect exact equality, not EXPECT_DOUBLE_EQ: this is the same ref twice
        EXPECT_EQ(ddc::get<Tag1>(val), ddc::get<Tag1>(chunk_cref(ix)));
        EXPECT_EQ(ddc::get<Tag2>(val), ddc::get<Tag2>(chunk_cref(ix)));
    }
}

TEST(VectorField1DTest, Access)
{
    Coord constexpr factor(1.391, 2.444);
    VectorFieldX<double> chunk(dom_x);
    for (auto&& ix : chunk.domain()) {
        Coord val = double(ix.uid()) * factor;
        ddcHelper::get<Tag1>(chunk)(ix) = ddc::get<Tag1>(val);
        ddcHelper::get<Tag2>(chunk)(ix) = ddc::get<Tag2>(val);
        // we expect exact equality, not EXPECT_DOUBLE_EQ: this is the same ref twice
        EXPECT_EQ(ddc::get<Tag1>(val), ddc::get<Tag1>(chunk(ix)));
        EXPECT_EQ(ddc::get<Tag2>(val), ddc::get<Tag2>(chunk(ix)));
    }
}

TEST(VectorField1DTest, SpanCview)
{
    Coord constexpr factor(1.391, 2.444);
    VectorFieldX<double> chunk(dom_x);
    for (auto&& ix : chunk.domain()) {
        Coord val = double(ix.uid()) * factor;
        ddcHelper::get<Tag1>(chunk)(ix) = ddc::get<Tag1>(val);
        ddcHelper::get<Tag2>(chunk)(ix) = ddc::get<Tag2>(val);
        // we expect exact equality, not EXPECT_DOUBLE_EQ: this is the same ref twice
        EXPECT_EQ(ddc::get<Tag1>(val), ddc::get<Tag1>(chunk.span_cview()(ix)));
        EXPECT_EQ(ddc::get<Tag2>(val), ddc::get<Tag2>(chunk.span_cview()(ix)));
    }
}

TEST(VectorField1DTest, ViewConst)
{
    Coord constexpr factor(1.391, 2.444);
    VectorFieldX<double> chunk(dom_x);
    VectorFieldX<double> const& chunk_cref = chunk;
    for (auto&& ix : chunk.domain()) {
        Coord val = double(ix.uid()) * factor;
        ddcHelper::get<Tag1>(chunk)(ix) = ddc::get<Tag1>(val);
        ddcHelper::get<Tag2>(chunk)(ix) = ddc::get<Tag2>(val);
        // we expect exact equality, not EXPECT_DOUBLE_EQ: this is the same ref twice
        EXPECT_EQ(ddc::get<Tag1>(val), ddc::get<Tag1>(chunk_cref.span_view()(ix)));
        EXPECT_EQ(ddc::get<Tag2>(val), ddc::get<Tag2>(chunk_cref.span_view()(ix)));
    }
}

TEST(VectorField1DTest, View)
{
    Coord constexpr factor(1.391, 2.444);
    VectorFieldX<double> chunk(dom_x);
    VectorFieldSpanX<double> chunk_span = chunk.span_view();
    for (auto&& ix : chunk.domain()) {
        Coord val = double(ix.uid()) * factor;
        ddcHelper::get<Tag1>(chunk_span)(ix) = ddc::get<Tag1>(val);
        ddcHelper::get<Tag2>(chunk_span)(ix) = ddc::get<Tag2>(val);
        // we expect exact equality, not EXPECT_DOUBLE_EQ: this is the same ref twice
        EXPECT_EQ(ddc::get<Tag1>(val), ddc::get<Tag1>(chunk(ix)));
        EXPECT_EQ(ddc::get<Tag2>(val), ddc::get<Tag2>(chunk(ix)));
    }
}

// \}
// Functions inherited from VectorFieldCommon (and free functions implemented for it) \{

// constructors are hidden

// assignment operators are hidden

// no slicing (operator[]) in 1D

// access (operator()) operators are hidden

// TODO: accessor

TEST(VectorField1DTest, Rank)
{
    VectorFieldX<double> chunk(dom_x);
    EXPECT_EQ(VectorFieldX<double>::rank(), 1);
}

// TODO: stride

// swap is hidden

TEST(VectorField1DTest, Domain)
{
    VectorFieldX<double> chunk(dom_x);
    EXPECT_EQ(dom_x, chunk.domain());
}

TEST(VectorField1DTest, DomainX)
{
    VectorFieldX<double> chunk(dom_x);
    EXPECT_EQ(dom_x, chunk.domain<DDimX>());
}

// TODO: data_handle()

// TODO: internal_mdspan

// TODO: allocation_mdspan

TEST(VectorField1DTest, GetDomainX)
{
    VectorFieldX<double> chunk(dom_x);
    EXPECT_EQ(dom_x, ddcHelper::get_domain<DDimX>(chunk));
}


TEST(VectorField1DTest, Deepcopy)
{
    VectorFieldX<double> chunk(dom_x);
    for (auto&& ix : chunk.domain()) {
        ddcHelper::get<Tag1>(chunk)(ix) = 1.001 * ix.uid();
        ddcHelper::get<Tag2>(chunk)(ix) = 0.0;
    }
    VectorFieldX<double> chunk2(chunk.domain());
    ddcHelper::deepcopy(chunk2, chunk);
    for (auto&& ix : chunk.domain()) {
        // we expect exact equality, not EXPECT_DOUBLE_EQ: these are copy
        EXPECT_EQ(ddc::get<Tag1>(chunk2(ix)), ddc::get<Tag1>(chunk(ix)));
        EXPECT_EQ(ddc::get<Tag2>(chunk2(ix)), ddc::get<Tag2>(chunk(ix)));
    }
}

// \}
// Functions implemented in VectorField 2D (and free functions specific to it) \{

// TODO: lots to do still!

TEST(VectorField2DTest, Access)
{
    VectorFieldXY<double> chunk(dom_x_y);
    for (auto&& ix : chunk.domain<DDimX>()) {
        for (auto&& iy : chunk.domain<DDimY>()) {
            ddcHelper::get<Tag1>(chunk)(ix, iy) = 1.357 * ix.uid();
            ddcHelper::get<Tag2>(chunk)(ix, iy) = 1.159 * iy.uid();
            // we expect exact equality, not EXPECT_DOUBLE_EQ: this is the same ref twice
            EXPECT_EQ(ddc::get<Tag1>(chunk(ix, iy)), ddc::get<Tag1>(chunk(ix, iy)));
            EXPECT_EQ(ddc::get<Tag2>(chunk(ix, iy)), ddc::get<Tag2>(chunk(ix, iy)));
        }
    }
}

TEST(VectorField2DTest, AccessReordered)
{
    VectorFieldXY<double> chunk(dom_x_y);
    for (auto&& ix : chunk.domain<DDimX>()) {
        for (auto&& iy : chunk.domain<DDimY>()) {
            ddcHelper::get<Tag1>(chunk)(ix, iy) = 1.455 * ix.uid();
            ddcHelper::get<Tag2>(chunk)(ix, iy) = 1.522 * iy.uid();
            // we expect exact equality, not EXPECT_DOUBLE_EQ: this is the same ref twice
            EXPECT_EQ(ddc::get<Tag1>(chunk(iy, ix)), ddc::get<Tag1>(chunk(ix, iy)));
            EXPECT_EQ(ddc::get<Tag2>(chunk(iy, ix)), ddc::get<Tag2>(chunk(ix, iy)));
        }
    }
}

TEST(VectorField2DTest, Cview)
{
    VectorFieldXY<double> chunk(dom_x_y);
    for (auto&& ix : chunk.domain<DDimX>()) {
        for (auto&& iy : chunk.domain<DDimY>()) {
            ddcHelper::get<Tag1>(chunk)(ix, iy) = 1. * ix.uid();
            ddcHelper::get<Tag2>(chunk)(ix, iy) = .001 * iy.uid();
        }
    }
    auto cview = chunk.span_cview();
    for (auto&& ix : chunk.domain<DDimX>()) {
        for (auto&& iy : chunk.domain<DDimY>()) {
            // we expect complete equality, not EXPECT_DOUBLE_EQ: these are copy
            EXPECT_EQ(ddc::get<Tag1>(cview(ix, iy)), ddc::get<Tag1>(chunk(ix, iy)));
            EXPECT_EQ(ddc::get<Tag2>(cview(ix, iy)), ddc::get<Tag2>(chunk(ix, iy)));
        }
    }
}

TEST(VectorField2DTest, SliceCoordX)
{
    DElemX constexpr slice_x_val = DElemX(lbound_x + 1);

    VectorFieldXY<double> chunk(dom_x_y);
    VectorFieldXY<double> const& chunk_cref = chunk;
    for (auto&& ix : chunk.domain<DDimX>()) {
        for (auto&& iy : chunk.domain<DDimY>()) {
            ddcHelper::get<Tag1>(chunk)(ix, iy) = 1. * ix.uid();
            ddcHelper::get<Tag2>(chunk)(ix, iy) = .001 * iy.uid();
        }
    }

    auto&& chunk_y = chunk_cref[slice_x_val];
    EXPECT_TRUE((std::is_same_v<
                 std::decay_t<decltype(chunk_y)>::layout_type,
                 std::experimental::layout_right>));
    EXPECT_EQ(
            ddcHelper::get<Tag1>(chunk_y).extent<DDimY>(),
            ddcHelper::get<Tag1>(chunk).extent<DDimY>());
    EXPECT_EQ(
            ddcHelper::get<Tag2>(chunk_y).extent<DDimY>(),
            ddcHelper::get<Tag2>(chunk).extent<DDimY>());
    for (auto&& ix : chunk_cref.domain<DDimY>()) {
        // we expect complete equality, not EXPECT_DOUBLE_EQ: these are copy
        EXPECT_EQ(ddc::get<Tag1>(chunk_y(ix)), ddc::get<Tag1>(chunk_cref(slice_x_val, ix)));
        EXPECT_EQ(ddc::get<Tag2>(chunk_y(ix)), ddc::get<Tag2>(chunk_cref(slice_x_val, ix)));
    }
}

TEST(VectorField2DTest, SliceCoordY)
{
    DElemY constexpr slice_y_val = DElemY(lbound_y + 1);

    VectorFieldXY<double> chunk(dom_x_y);
    VectorFieldXY<double> const& chunk_cref = chunk;
    for (auto&& ix : chunk.domain<DDimX>()) {
        for (auto&& iy : chunk.domain<DDimY>()) {
            ddcHelper::get<Tag1>(chunk)(ix, iy) = 1. * ix.uid();
            ddcHelper::get<Tag2>(chunk)(ix, iy) = .001 * iy.uid();
        }
    }

    auto&& chunk_x = chunk_cref[slice_y_val];
    EXPECT_TRUE((std::is_same_v<
                 std::decay_t<decltype(chunk_x)>::layout_type,
                 std::experimental::layout_stride>));
    EXPECT_EQ(
            ddcHelper::get<Tag1>(chunk_x).extent<DDimX>(),
            ddcHelper::get<Tag1>(chunk).extent<DDimX>());
    EXPECT_EQ(
            ddcHelper::get<Tag2>(chunk_x).extent<DDimX>(),
            ddcHelper::get<Tag2>(chunk).extent<DDimX>());
    for (auto&& ix : chunk_cref.domain<DDimX>()) {
        // we expect complete equality, not EXPECT_DOUBLE_EQ: these are copy
        EXPECT_EQ(ddc::get<Tag1>(chunk_x(ix)), ddc::get<Tag1>(chunk_cref(ix, slice_y_val)));
        EXPECT_EQ(ddc::get<Tag2>(chunk_x(ix)), ddc::get<Tag2>(chunk_cref(ix, slice_y_val)));
    }
}

TEST(VectorField2DTest, SliceDomainX)
{
    DDomX constexpr subdomain_x = DDomX(DElemX(lbound_x + 1), DVectX(nelems_x - 2));

    VectorFieldXY<double> chunk(dom_x_y);
    VectorFieldXY<double> const& chunk_cref = chunk;
    for (auto&& ix : chunk.domain<DDimX>()) {
        for (auto&& iy : chunk.domain<DDimY>()) {
            ddcHelper::get<Tag1>(chunk)(ix, iy) = 1. * ix.uid();
            ddcHelper::get<Tag2>(chunk)(ix, iy) = .001 * iy.uid();
        }
    }

    auto&& subchunk_x = chunk_cref[subdomain_x];
    EXPECT_TRUE((std::is_same_v<
                 std::decay_t<decltype(subchunk_x)>::layout_type,
                 std::experimental::layout_right>));

    EXPECT_EQ(ddcHelper::get<Tag1>(subchunk_x).extent<DDimX>(), subdomain_x.size());
    EXPECT_EQ(ddcHelper::get<Tag2>(subchunk_x).extent<DDimY>(), chunk.domain<DDimY>().size());
    EXPECT_EQ(ddcHelper::get<Tag1>(subchunk_x).extent<DDimX>(), subdomain_x.size());
    EXPECT_EQ(ddcHelper::get<Tag2>(subchunk_x).extent<DDimY>(), chunk.domain<DDimY>().size());
    for (auto&& ix : subchunk_x.domain<DDimX>()) {
        for (auto&& iy : subchunk_x.domain<DDimY>()) {
            // we expect complete equality, not EXPECT_DOUBLE_EQ: these are copy
            EXPECT_EQ(ddc::get<Tag1>(subchunk_x(ix, iy)), ddc::get<Tag1>(chunk_cref(ix, iy)));
            EXPECT_EQ(ddc::get<Tag2>(subchunk_x(ix, iy)), ddc::get<Tag2>(chunk_cref(ix, iy)));
        }
    }
}

TEST(VectorField2DTest, SliceDomainXTooearly)
{
    [[maybe_unused]] DDomX constexpr subdomain_x = DDomX(DElemX(lbound_x - 1), nelems_x);

    VectorFieldXY<double> chunk(dom_x_y);
#ifndef NDEBUG // The assertion is only checked if NDEBUG isn't defined
    // the error message is checked with clang & gcc only
    EXPECT_DEATH(
            chunk[subdomain_x],
            R"rgx([Aa]ssert.*uid<ODDims>\(m_element_begin\).*uid<ODDims>\(odomain\.m_element_begin\))rgx");
#endif
}

TEST(VectorField2DTest, SliceDomainXToolate)
{
    [[maybe_unused]] DDomX constexpr subdomain_x = DDomX(lbound_x, DVectX(nelems_x + 1));

    VectorFieldXY<double> chunk(dom_x_y);
#ifndef NDEBUG // The assertion is only checked if NDEBUG isn't defined
    // the error message is checked with clang & gcc only
    EXPECT_DEATH(
            chunk[subdomain_x],
            R"rgx([Aa]ssert.*uid<ODDims>\(m_element_end\).*uid<ODDims>\(odomain\.m_element_end\).*)rgx");
#endif
}

TEST(VectorField2DTest, SliceDomainY)
{
    DDomY constexpr subdomain_y = DDomY(DElemY(lbound_y + 1), DVectY(nelems_y - 2));

    VectorFieldXY<double> chunk(dom_x_y);
    VectorFieldXY<double> const& chunk_cref = chunk;
    for (auto&& ix : chunk.domain<DDimX>()) {
        for (auto&& iy : chunk.domain<DDimY>()) {
            ddcHelper::get<Tag1>(chunk)(ix, iy) = 1. * ix.uid();
            ddcHelper::get<Tag2>(chunk)(ix, iy) = .001 * iy.uid();
        }
    }
    auto&& subchunk_y = chunk_cref[subdomain_y];
    EXPECT_TRUE((std::is_same_v<
                 std::decay_t<decltype(subchunk_y)>::layout_type,
                 std::experimental::layout_stride>));

    EXPECT_EQ(ddcHelper::get<Tag1>(subchunk_y).extent<DDimX>(), chunk.domain<DDimX>().size());
    EXPECT_EQ(ddcHelper::get<Tag1>(subchunk_y).extent<DDimY>(), subdomain_y.size());
    EXPECT_EQ(ddcHelper::get<Tag2>(subchunk_y).extent<DDimX>(), chunk.domain<DDimX>().size());
    EXPECT_EQ(ddcHelper::get<Tag2>(subchunk_y).extent<DDimY>(), subdomain_y.size());
    for (auto&& ix : subchunk_y.domain<DDimX>()) {
        for (auto&& iy : subchunk_y.domain<DDimY>()) {
            // we expect complete equality, not EXPECT_DOUBLE_EQ: these are copy
            EXPECT_EQ(ddc::get<Tag1>(subchunk_y(ix, iy)), ddc::get<Tag1>(chunk_cref(ix, iy)));
            EXPECT_EQ(ddc::get<Tag2>(subchunk_y(ix, iy)), ddc::get<Tag2>(chunk_cref(ix, iy)));
        }
    }
}

TEST(VectorField2DTest, Deepcopy)
{
    VectorFieldXY<double> chunk(dom_x_y);
    for (auto&& ix : chunk.domain<DDimX>()) {
        for (auto&& iy : chunk.domain<DDimY>()) {
            ddcHelper::get<Tag1>(chunk)(ix, iy) = 1.739 * ix.uid();
            ddcHelper::get<Tag2>(chunk)(ix, iy) = 1.412 * iy.uid();
        }
    }
    VectorFieldXY<double> chunk2(chunk.domain());
    ddcHelper::deepcopy(chunk2, chunk);
    for (auto&& ix : chunk.domain<DDimX>()) {
        for (auto&& iy : chunk.domain<DDimY>()) {
            // we expect complete equality, not EXPECT_DOUBLE_EQ: these are copy
            EXPECT_EQ(ddc::get<Tag1>(chunk2(ix, iy)), ddc::get<Tag1>(chunk(ix, iy)));
            EXPECT_EQ(ddc::get<Tag2>(chunk2(ix, iy)), ddc::get<Tag2>(chunk(ix, iy)));
        }
    }
}

TEST(VectorField2DTest, DeepcopyReordered)
{
    VectorFieldXY<double> chunk(dom_x_y);
    for (auto&& ix : chunk.domain<DDimX>()) {
        for (auto&& iy : chunk.domain<DDimY>()) {
            ddcHelper::get<Tag1>(chunk)(ix, iy) = 1.739 * ix.uid();
            ddcHelper::get<Tag2>(chunk)(ix, iy) = 1.412 * iy.uid();
        }
    }
    VectorFieldYX<double> chunk2(ddc::select<DDimY, DDimX>(chunk.domain()));
    VectorFieldSpan<double, DDomXY, Direction, std::experimental::layout_left> chunk2_view(
            chunk.domain(),
            chunk2.get<Tag1>().data_handle(),
            chunk2.get<Tag2>().data_handle());
    ddcHelper::deepcopy(chunk2_view, chunk);
    for (auto&& ix : chunk.domain<DDimX>()) {
        for (auto&& iy : chunk.domain<DDimY>()) {
            // we expect complete equality, not EXPECT_DOUBLE_EQ: these are copy
            EXPECT_EQ(ddc::get<Tag1>(chunk2(ix, iy)), ddc::get<Tag1>(chunk(ix, iy)));
            EXPECT_EQ(ddc::get<Tag2>(chunk2(ix, iy)), ddc::get<Tag2>(chunk(ix, iy)));
            EXPECT_EQ(ddc::get<Tag1>(chunk2(ix, iy)), ddc::get<Tag1>(chunk(iy, ix)));
            EXPECT_EQ(ddc::get<Tag2>(chunk2(ix, iy)), ddc::get<Tag2>(chunk(iy, ix)));
        }
    }
}
