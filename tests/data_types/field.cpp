// SPDX-License-Identifier: MIT
#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include "ddc_aliases.hpp"
#include "directional_tag.hpp"
#include "vector_field.hpp"
#include "vector_field_span.hpp"

namespace {

class Tag1
{
};
class Tag2
{
};

using Direction = NDTag<Tag1, Tag2>;

using Coord2D = Coord<Tag1, Tag2>;

using Idx0D = Idx<>;
using IdxStept0D = IdxStep<>;
using IdxRange0D = IdxRange<>;

template <class Datatype>
using VectorField0D = VectorField<Datatype, IdxRange0D, Direction>;
template <class Datatype>
using VectorFieldField0D = VectorFieldSpan<Datatype, IdxRange0D, Direction>;

struct GridX
{
};
using IdxX = Idx<GridX>;
using IdxSteptX = IdxStep<GridX>;
using IdxRangeX = IdxRange<GridX>;

template <class Datatype>
using VectorFieldX = VectorField<Datatype, IdxRangeX, Direction>;
template <class Datatype>
using VectorFieldFieldX = VectorFieldSpan<Datatype, IdxRangeX, Direction>;


struct GridY
{
};
using IdxY = Idx<GridY>;
using IdxSteptY = IdxStep<GridY>;
using IdxRangeY = IdxRange<GridY>;

template <class Datatype>
using VectorFieldY = VectorField<Datatype, IdxRangeY, Direction>;
template <class Datatype>
using VectorFieldFieldY = VectorFieldSpan<Datatype, IdxRangeY, Direction>;


struct GridZ
{
};
using IdxZ = Idx<GridZ>;
using IdxSteptZ = IdxStep<GridZ>;
using IdxRangeZ = IdxRange<GridZ>;


using IdxXY = Idx<GridX, GridY>;
using IdxSteptXY = IdxStep<GridX, GridY>;
using IdxRangeXY = IdxRange<GridX, GridY>;

template <class Datatype>
using VectorFieldXY = VectorField<Datatype, IdxRangeXY, Direction>;


using IdxYX = Idx<GridY, GridX>;
using IdxSteptYX = IdxStep<GridY, GridX>;
using IdxRangeYX = IdxRange<GridY, GridX>;

template <class Datatype>
using VectorFieldYX = VectorField<Datatype, IdxRangeYX, Direction>;
template <class Datatype>
using VectorFieldFieldYX = VectorFieldSpan<Datatype, IdxRangeYX, Direction>;


static Idx0D constexpr lbound_0d {};
static IdxStept0D constexpr nelems_0d {};
static IdxRange0D constexpr dom_0d(lbound_0d, nelems_0d);

static IdxX constexpr lbound_x(50);
static IdxSteptX constexpr nelems_x(3);
static IdxRangeX constexpr dom_x(lbound_x, nelems_x);

static IdxY constexpr lbound_y(4);
static IdxSteptY constexpr nelems_y(12);

static IdxXY constexpr lbound_x_y {lbound_x, lbound_y};
static IdxSteptXY constexpr nelems_x_y(nelems_x, nelems_y);
static IdxRangeXY constexpr dom_x_y(lbound_x_y, nelems_x_y);

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
    Coord2D constexpr factor(1.391, 2.444);
    VectorField0D<double> chunk(dom_0d);
    ddcHelper::get<Tag1>(chunk)() = ddc::get<Tag1>(factor);
    ddcHelper::get<Tag2>(chunk)() = ddc::get<Tag2>(factor);

    VectorField0D<double> chunk2(std::move(chunk));
    EXPECT_EQ(get_idx_range(chunk2), dom_0d);
    EXPECT_DOUBLE_EQ(ddc::get<Tag1>(factor), ddc::get<Tag1>(chunk2()));
    EXPECT_DOUBLE_EQ(ddc::get<Tag2>(factor), ddc::get<Tag2>(chunk2()));
}

TEST(VectorField1DTest, MoveConstructor)
{
    Coord2D constexpr factor(1.391, 2.444);
    VectorFieldX<double> chunk(dom_x);
    for (auto&& ix : get_idx_range(chunk)) {
        Coord2D val = double(ix.uid()) * factor;
        ddcHelper::get<Tag1>(chunk)(ix) = ddc::get<Tag1>(val);
        ddcHelper::get<Tag2>(chunk)(ix) = ddc::get<Tag2>(val);
    }

    VectorFieldX<double> chunk2(std::move(chunk));
    EXPECT_EQ(get_idx_range(chunk2), dom_x);
    for (auto&& ix : get_idx_range(chunk2)) {
        // we expect exact equality, not EXPECT_DOUBLE_EQ: this is the same ref twice
        Coord2D val = double(ix.uid()) * factor;
        EXPECT_EQ(ddc::get<Tag1>(val), ddc::get<Tag1>(chunk2(ix)));
        EXPECT_EQ(ddc::get<Tag2>(val), ddc::get<Tag2>(chunk2(ix)));
    }
}

TEST(VectorField0DTest, MoveAssignment)
{
    Coord2D constexpr factor(1.391, 2.444);
    VectorField0D<double> chunk(dom_0d);
    ddcHelper::get<Tag1>(chunk)() = ddc::get<Tag1>(factor);
    ddcHelper::get<Tag2>(chunk)() = ddc::get<Tag2>(factor);

    VectorField0D<double> chunk2(IdxRange0D(lbound_0d, IdxStept0D()));
    chunk2 = std::move(chunk);
    EXPECT_EQ(get_idx_range(chunk2), dom_0d);
    EXPECT_DOUBLE_EQ(ddc::get<Tag1>(factor), ddc::get<Tag1>(chunk2()));
    EXPECT_DOUBLE_EQ(ddc::get<Tag2>(factor), ddc::get<Tag2>(chunk2()));
}

TEST(VectorField1DTest, MoveAssignment)
{
    Coord2D constexpr factor(1.391, 2.444);
    VectorFieldX<double> chunk(dom_x);
    for (auto&& ix : get_idx_range(chunk)) {
        Coord2D val = double(ix.uid()) * factor;
        ddcHelper::get<Tag1>(chunk)(ix) = ddc::get<Tag1>(val);
        ddcHelper::get<Tag2>(chunk)(ix) = ddc::get<Tag2>(val);
    }

    VectorFieldX<double> chunk2(IdxRangeX(lbound_x, IdxSteptX(0)));
    chunk2 = std::move(chunk);
    EXPECT_EQ(get_idx_range(chunk2), dom_x);
    for (auto&& ix : get_idx_range(chunk2)) {
        // we expect exact equality, not EXPECT_DOUBLE_EQ: this is the same ref twice
        Coord2D val = double(ix.uid()) * factor;
        EXPECT_EQ(ddc::get<Tag1>(val), ddc::get<Tag1>(chunk2(ix)));
        EXPECT_EQ(ddc::get<Tag2>(val), ddc::get<Tag2>(chunk2(ix)));
    }
}

TEST(VectorField0DTest, Swap)
{
    Coord2D constexpr factor(1.391, 2.444);
    VectorField0D<double> chunk(dom_0d);
    ddcHelper::get<Tag1>(chunk)() = ddc::get<Tag1>(factor);
    ddcHelper::get<Tag2>(chunk)() = ddc::get<Tag2>(factor);

    IdxRange0D empty_idx_range(lbound_0d, IdxStept0D());
    VectorField0D<double> chunk2(empty_idx_range);

    std::swap(chunk2, chunk);
    EXPECT_EQ(get_idx_range(chunk), empty_idx_range);
    EXPECT_EQ(get_idx_range(chunk2), dom_0d);
    EXPECT_DOUBLE_EQ(ddc::get<Tag1>(factor), ddc::get<Tag1>(chunk2()));
    EXPECT_DOUBLE_EQ(ddc::get<Tag2>(factor), ddc::get<Tag2>(chunk2()));
}

TEST(VectorField1DTest, Swap)
{
    Coord2D constexpr factor(1.391, 2.444);
    VectorFieldX<double> chunk(dom_x);
    for (auto&& ix : get_idx_range(chunk)) {
        Coord2D val = double(ix.uid()) * factor;
        ddcHelper::get<Tag1>(chunk)(ix) = ddc::get<Tag1>(val);
        ddcHelper::get<Tag2>(chunk)(ix) = ddc::get<Tag2>(val);
    }

    IdxRangeX empty_idx_range(lbound_x, IdxSteptX(0));
    VectorFieldX<double> chunk2(empty_idx_range);

    std::swap(chunk2, chunk);
    EXPECT_EQ(get_idx_range(chunk), empty_idx_range);
    EXPECT_EQ(get_idx_range(chunk2), dom_x);
    for (auto&& ix : get_idx_range(chunk2)) {
        // we expect exact equality, not EXPECT_DOUBLE_EQ: this is the same ref twice
        Coord2D val = double(ix.uid()) * factor;
        EXPECT_EQ(ddc::get<Tag1>(val), ddc::get<Tag1>(chunk2(ix)));
        EXPECT_EQ(ddc::get<Tag2>(val), ddc::get<Tag2>(chunk2(ix)));
    }
}

// no dim subset access in 1D

TEST(VectorField1DTest, AccessConst)
{
    Coord2D constexpr factor(1.391, 2.444);
    VectorFieldX<double> chunk(dom_x);
    VectorFieldX<double> const& chunk_cref = chunk;
    for (auto&& ix : get_idx_range(chunk)) {
        Coord2D val = double(ix.uid()) * factor;
        ddcHelper::get<Tag1>(chunk)(ix) = ddc::get<Tag1>(val);
        ddcHelper::get<Tag2>(chunk)(ix) = ddc::get<Tag2>(val);
        // we expect exact equality, not EXPECT_DOUBLE_EQ: this is the same ref twice
        EXPECT_EQ(ddc::get<Tag1>(val), ddc::get<Tag1>(chunk_cref(ix)));
        EXPECT_EQ(ddc::get<Tag2>(val), ddc::get<Tag2>(chunk_cref(ix)));
    }
}

TEST(VectorField1DTest, Access)
{
    Coord2D constexpr factor(1.391, 2.444);
    VectorFieldX<double> chunk(dom_x);
    for (auto&& ix : get_idx_range(chunk)) {
        Coord2D val = double(ix.uid()) * factor;
        ddcHelper::get<Tag1>(chunk)(ix) = ddc::get<Tag1>(val);
        ddcHelper::get<Tag2>(chunk)(ix) = ddc::get<Tag2>(val);
        // we expect exact equality, not EXPECT_DOUBLE_EQ: this is the same ref twice
        EXPECT_EQ(ddc::get<Tag1>(val), ddc::get<Tag1>(chunk(ix)));
        EXPECT_EQ(ddc::get<Tag2>(val), ddc::get<Tag2>(chunk(ix)));
    }
}

TEST(VectorField1DTest, SpanCview)
{
    Coord2D constexpr factor(1.391, 2.444);
    VectorFieldX<double> chunk(dom_x);
    for (auto&& ix : get_idx_range(chunk)) {
        Coord2D val = double(ix.uid()) * factor;
        ddcHelper::get<Tag1>(chunk)(ix) = ddc::get<Tag1>(val);
        ddcHelper::get<Tag2>(chunk)(ix) = ddc::get<Tag2>(val);
        // we expect exact equality, not EXPECT_DOUBLE_EQ: this is the same ref twice
        EXPECT_EQ(ddc::get<Tag1>(val), ddc::get<Tag1>(get_const_field(chunk)(ix)));
        EXPECT_EQ(ddc::get<Tag2>(val), ddc::get<Tag2>(get_const_field(chunk)(ix)));
    }
}

TEST(VectorField1DTest, ViewConst)
{
    Coord2D constexpr factor(1.391, 2.444);
    VectorFieldX<double> chunk(dom_x);
    VectorFieldX<double> const& chunk_cref = chunk;
    for (auto&& ix : get_idx_range(chunk)) {
        Coord2D val = double(ix.uid()) * factor;
        ddcHelper::get<Tag1>(chunk)(ix) = ddc::get<Tag1>(val);
        ddcHelper::get<Tag2>(chunk)(ix) = ddc::get<Tag2>(val);
        // we expect exact equality, not EXPECT_DOUBLE_EQ: this is the same ref twice
        EXPECT_EQ(ddc::get<Tag1>(val), ddc::get<Tag1>(get_field(chunk_cref)(ix)));
        EXPECT_EQ(ddc::get<Tag2>(val), ddc::get<Tag2>(get_field(chunk_cref)(ix)));
    }
}

TEST(VectorField1DTest, View)
{
    Coord2D constexpr factor(1.391, 2.444);
    VectorFieldX<double> chunk(dom_x);
    VectorFieldFieldX<double> chunk_span = get_field(chunk);
    for (auto&& ix : get_idx_range(chunk)) {
        Coord2D val = double(ix.uid()) * factor;
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
    EXPECT_EQ(dom_x, get_idx_range(chunk));
}

TEST(VectorField1DTest, DomainX)
{
    VectorFieldX<double> chunk(dom_x);
    EXPECT_EQ(dom_x, get_idx_range<GridX>(chunk));
}

// TODO: data_handle()

// TODO: internal_mdspan

// TODO: allocation_mdspan

TEST(VectorField1DTest, GetDomainX)
{
    VectorFieldX<double> chunk(dom_x);
    EXPECT_EQ(dom_x, ddcHelper::get_domain<GridX>(chunk));
}


TEST(VectorField1DTest, Deepcopy)
{
    VectorFieldX<double> chunk(dom_x);
    for (auto&& ix : get_idx_range(chunk)) {
        ddcHelper::get<Tag1>(chunk)(ix) = 1.001 * ix.uid();
        ddcHelper::get<Tag2>(chunk)(ix) = 0.0;
    }
    VectorFieldX<double> chunk2(get_idx_range(chunk));
    ddcHelper::deepcopy(chunk2, chunk);
    for (auto&& ix : get_idx_range(chunk)) {
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
    for (auto&& ix : get_idx_range<GridX>(chunk)) {
        for (auto&& iy : get_idx_range<GridY>(chunk)) {
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
    for (auto&& ix : get_idx_range<GridX>(chunk)) {
        for (auto&& iy : get_idx_range<GridY>(chunk)) {
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
    for (auto&& ix : get_idx_range<GridX>(chunk)) {
        for (auto&& iy : get_idx_range<GridY>(chunk)) {
            ddcHelper::get<Tag1>(chunk)(ix, iy) = 1. * ix.uid();
            ddcHelper::get<Tag2>(chunk)(ix, iy) = .001 * iy.uid();
        }
    }
    auto cview = get_const_field(chunk);
    for (auto&& ix : get_idx_range<GridX>(chunk)) {
        for (auto&& iy : get_idx_range<GridY>(chunk)) {
            // we expect complete equality, not EXPECT_DOUBLE_EQ: these are copy
            EXPECT_EQ(ddc::get<Tag1>(cview(ix, iy)), ddc::get<Tag1>(chunk(ix, iy)));
            EXPECT_EQ(ddc::get<Tag2>(cview(ix, iy)), ddc::get<Tag2>(chunk(ix, iy)));
        }
    }
}

TEST(VectorField2DTest, SliceCoordX)
{
    IdxX constexpr slice_x_val = IdxX(lbound_x + 1);

    VectorFieldXY<double> chunk(dom_x_y);
    VectorFieldXY<double> const& chunk_cref = chunk;
    for (auto&& ix : get_idx_range<GridX>(chunk)) {
        for (auto&& iy : get_idx_range<GridY>(chunk)) {
            ddcHelper::get<Tag1>(chunk)(ix, iy) = 1. * ix.uid();
            ddcHelper::get<Tag2>(chunk)(ix, iy) = .001 * iy.uid();
        }
    }

    auto&& chunk_y = chunk_cref[slice_x_val];
    EXPECT_TRUE((std::is_same_v<
                 std::decay_t<decltype(chunk_y)>::layout_type,
                 std::experimental::layout_right>));
    EXPECT_EQ(
            ddcHelper::get<Tag1>(chunk_y).extent<GridY>(),
            ddcHelper::get<Tag1>(chunk).extent<GridY>());
    EXPECT_EQ(
            ddcHelper::get<Tag2>(chunk_y).extent<GridY>(),
            ddcHelper::get<Tag2>(chunk).extent<GridY>());
    for (auto&& ix : get_idx_range<GridY>(chunk_cref)) {
        // we expect complete equality, not EXPECT_DOUBLE_EQ: these are copy
        EXPECT_EQ(ddc::get<Tag1>(chunk_y(ix)), ddc::get<Tag1>(chunk_cref(slice_x_val, ix)));
        EXPECT_EQ(ddc::get<Tag2>(chunk_y(ix)), ddc::get<Tag2>(chunk_cref(slice_x_val, ix)));
    }
}

TEST(VectorField2DTest, SliceCoordY)
{
    IdxY constexpr slice_y_val = IdxY(lbound_y + 1);

    VectorFieldXY<double> chunk(dom_x_y);
    VectorFieldXY<double> const& chunk_cref = chunk;
    for (auto&& ix : get_idx_range<GridX>(chunk)) {
        for (auto&& iy : get_idx_range<GridY>(chunk)) {
            ddcHelper::get<Tag1>(chunk)(ix, iy) = 1. * ix.uid();
            ddcHelper::get<Tag2>(chunk)(ix, iy) = .001 * iy.uid();
        }
    }

    auto&& chunk_x = chunk_cref[slice_y_val];
    EXPECT_TRUE((std::is_same_v<
                 std::decay_t<decltype(chunk_x)>::layout_type,
                 std::experimental::layout_stride>));
    EXPECT_EQ(
            ddcHelper::get<Tag1>(chunk_x).extent<GridX>(),
            ddcHelper::get<Tag1>(chunk).extent<GridX>());
    EXPECT_EQ(
            ddcHelper::get<Tag2>(chunk_x).extent<GridX>(),
            ddcHelper::get<Tag2>(chunk).extent<GridX>());
    for (auto&& ix : get_idx_range<GridX>(chunk_cref)) {
        // we expect complete equality, not EXPECT_DOUBLE_EQ: these are copy
        EXPECT_EQ(ddc::get<Tag1>(chunk_x(ix)), ddc::get<Tag1>(chunk_cref(ix, slice_y_val)));
        EXPECT_EQ(ddc::get<Tag2>(chunk_x(ix)), ddc::get<Tag2>(chunk_cref(ix, slice_y_val)));
    }
}

TEST(VectorField2DTest, SliceDomainX)
{
    IdxRangeX constexpr subidx_range_x = IdxRangeX(IdxX(lbound_x + 1), IdxSteptX(nelems_x - 2));

    VectorFieldXY<double> chunk(dom_x_y);
    VectorFieldXY<double> const& chunk_cref = chunk;
    for (auto&& ix : get_idx_range<GridX>(chunk)) {
        for (auto&& iy : get_idx_range<GridY>(chunk)) {
            ddcHelper::get<Tag1>(chunk)(ix, iy) = 1. * ix.uid();
            ddcHelper::get<Tag2>(chunk)(ix, iy) = .001 * iy.uid();
        }
    }

    auto&& subchunk_x = chunk_cref[subidx_range_x];
    EXPECT_TRUE((std::is_same_v<
                 std::decay_t<decltype(subchunk_x)>::layout_type,
                 std::experimental::layout_right>));

    EXPECT_EQ(ddcHelper::get<Tag1>(subchunk_x).extent<GridX>(), subidx_range_x.size());
    EXPECT_EQ(ddcHelper::get<Tag2>(subchunk_x).extent<GridY>(), get_idx_range<GridY>(chunk).size());
    EXPECT_EQ(ddcHelper::get<Tag1>(subchunk_x).extent<GridX>(), subidx_range_x.size());
    EXPECT_EQ(ddcHelper::get<Tag2>(subchunk_x).extent<GridY>(), get_idx_range<GridY>(chunk).size());
    for (auto&& ix : get_idx_range<GridX>(subchunk_x)) {
        for (auto&& iy : get_idx_range<GridY>(subchunk_x)) {
            // we expect complete equality, not EXPECT_DOUBLE_EQ: these are copy
            EXPECT_EQ(ddc::get<Tag1>(subchunk_x(ix, iy)), ddc::get<Tag1>(chunk_cref(ix, iy)));
            EXPECT_EQ(ddc::get<Tag2>(subchunk_x(ix, iy)), ddc::get<Tag2>(chunk_cref(ix, iy)));
        }
    }
}

TEST(VectorField2DTest, SliceDomainXTooearly)
{
    [[maybe_unused]] IdxRangeX constexpr subidx_range_x = IdxRangeX(IdxX(lbound_x - 1), nelems_x);

    VectorFieldXY<double> chunk(dom_x_y);
#ifndef NDEBUG // The assertion is only checked if NDEBUG isn't defined
    // the error message is checked with clang & gcc only
    EXPECT_DEATH(
            chunk[subidx_range_x],
            R"rgx([Aa]ssert.*uid<ODDims>\(m_element_begin\).*uid<ODDims>\(odomain\.m_element_begin\))rgx");
#endif
}

TEST(VectorField2DTest, SliceDomainXToolate)
{
    [[maybe_unused]] IdxRangeX constexpr subidx_range_x
            = IdxRangeX(lbound_x, IdxSteptX(nelems_x + 1));

    VectorFieldXY<double> chunk(dom_x_y);
#ifndef NDEBUG // The assertion is only checked if NDEBUG isn't defined
    // the error message is checked with clang & gcc only
    EXPECT_DEATH(
            chunk[subidx_range_x],
            R"rgx([Aa]ssert.*uid<ODDims>\(m_element_end\).*uid<ODDims>\(odomain\.m_element_end\).*)rgx");
#endif
}

TEST(VectorField2DTest, SliceDomainY)
{
    IdxRangeY constexpr subidx_range_y = IdxRangeY(IdxY(lbound_y + 1), IdxSteptY(nelems_y - 2));

    VectorFieldXY<double> chunk(dom_x_y);
    VectorFieldXY<double> const& chunk_cref = chunk;
    for (auto&& ix : get_idx_range<GridX>(chunk)) {
        for (auto&& iy : get_idx_range<GridY>(chunk)) {
            ddcHelper::get<Tag1>(chunk)(ix, iy) = 1. * ix.uid();
            ddcHelper::get<Tag2>(chunk)(ix, iy) = .001 * iy.uid();
        }
    }
    auto&& subchunk_y = chunk_cref[subidx_range_y];
    EXPECT_TRUE((std::is_same_v<
                 std::decay_t<decltype(subchunk_y)>::layout_type,
                 std::experimental::layout_stride>));

    EXPECT_EQ(ddcHelper::get<Tag1>(subchunk_y).extent<GridX>(), get_idx_range<GridX>(chunk).size());
    EXPECT_EQ(ddcHelper::get<Tag1>(subchunk_y).extent<GridY>(), subidx_range_y.size());
    EXPECT_EQ(ddcHelper::get<Tag2>(subchunk_y).extent<GridX>(), get_idx_range<GridX>(chunk).size());
    EXPECT_EQ(ddcHelper::get<Tag2>(subchunk_y).extent<GridY>(), subidx_range_y.size());
    for (auto&& ix : get_idx_range<GridX>(subchunk_y)) {
        for (auto&& iy : get_idx_range<GridY>(subchunk_y)) {
            // we expect complete equality, not EXPECT_DOUBLE_EQ: these are copy
            EXPECT_EQ(ddc::get<Tag1>(subchunk_y(ix, iy)), ddc::get<Tag1>(chunk_cref(ix, iy)));
            EXPECT_EQ(ddc::get<Tag2>(subchunk_y(ix, iy)), ddc::get<Tag2>(chunk_cref(ix, iy)));
        }
    }
}

TEST(VectorField2DTest, Deepcopy)
{
    VectorFieldXY<double> chunk(dom_x_y);
    for (auto&& ix : get_idx_range<GridX>(chunk)) {
        for (auto&& iy : get_idx_range<GridY>(chunk)) {
            ddcHelper::get<Tag1>(chunk)(ix, iy) = 1.739 * ix.uid();
            ddcHelper::get<Tag2>(chunk)(ix, iy) = 1.412 * iy.uid();
        }
    }
    VectorFieldXY<double> chunk2(get_idx_range(chunk));
    ddcHelper::deepcopy(chunk2, chunk);
    for (auto&& ix : get_idx_range<GridX>(chunk)) {
        for (auto&& iy : get_idx_range<GridY>(chunk)) {
            // we expect complete equality, not EXPECT_DOUBLE_EQ: these are copy
            EXPECT_EQ(ddc::get<Tag1>(chunk2(ix, iy)), ddc::get<Tag1>(chunk(ix, iy)));
            EXPECT_EQ(ddc::get<Tag2>(chunk2(ix, iy)), ddc::get<Tag2>(chunk(ix, iy)));
        }
    }
}

TEST(VectorField2DTest, DeepcopyReordered)
{
    VectorFieldXY<double> chunk(dom_x_y);
    for (auto&& ix : get_idx_range<GridX>(chunk)) {
        for (auto&& iy : get_idx_range<GridY>(chunk)) {
            ddcHelper::get<Tag1>(chunk)(ix, iy) = 1.739 * ix.uid();
            ddcHelper::get<Tag2>(chunk)(ix, iy) = 1.412 * iy.uid();
        }
    }
    VectorFieldYX<double> chunk2(ddc::select<GridY, GridX>(get_idx_range(chunk)));
    VectorFieldSpan<double, IdxRangeXY, Direction, std::experimental::layout_left> chunk2_view(
            get_idx_range(chunk),
            chunk2.get<Tag1>().data_handle(),
            chunk2.get<Tag2>().data_handle());
    ddcHelper::deepcopy(chunk2_view, chunk);
    for (auto&& ix : get_idx_range<GridX>(chunk)) {
        for (auto&& iy : get_idx_range<GridY>(chunk)) {
            // we expect complete equality, not EXPECT_DOUBLE_EQ: these are copy
            EXPECT_EQ(ddc::get<Tag1>(chunk2(ix, iy)), ddc::get<Tag1>(chunk(ix, iy)));
            EXPECT_EQ(ddc::get<Tag2>(chunk2(ix, iy)), ddc::get<Tag2>(chunk(ix, iy)));
            EXPECT_EQ(ddc::get<Tag1>(chunk2(ix, iy)), ddc::get<Tag1>(chunk(iy, ix)));
            EXPECT_EQ(ddc::get<Tag2>(chunk2(ix, iy)), ddc::get<Tag2>(chunk(iy, ix)));
        }
    }
}
