// SPDX-License-Identifier: MIT
#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "directional_tag.hpp"
#include "vector_field.hpp"
#include "vector_field_mem.hpp"

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
using IdxStep0D = IdxStep<>;
using IdxRange0D = IdxRange<>;

using DVectorFieldMem0D = VectorFieldMem<double, IdxRange0D, Direction>;

struct GridX
{
};
using IdxX = Idx<GridX>;
using IdxStepX = IdxStep<GridX>;
using IdxRangeX = IdxRange<GridX>;

using DVectorFieldMemX = VectorFieldMem<double, IdxRangeX, Direction>;
using DVectorFieldX = VectorField<double, IdxRangeX, Direction>;
using DVectorConstFieldX = VectorField<double const, IdxRangeX, Direction>;
using DVectorConstFieldSliceX
        = VectorField<double const, IdxRangeX, Direction, std::experimental::layout_stride>;


struct GridY
{
};
using IdxY = Idx<GridY>;
using IdxStepY = IdxStep<GridY>;
using IdxRangeY = IdxRange<GridY>;

using DVectorFieldMemY = VectorFieldMem<double, IdxRangeY, Direction>;
using DVectorFieldY = VectorField<double, IdxRangeY, Direction>;
using DVectorConstFieldY = VectorField<double const, IdxRangeY, Direction>;


struct GridZ
{
};
using IdxZ = Idx<GridZ>;
using IdxStepZ = IdxStep<GridZ>;
using IdxRangeZ = IdxRange<GridZ>;


using IdxXY = Idx<GridX, GridY>;
using IdxStepXY = IdxStep<GridX, GridY>;
using IdxRangeXY = IdxRange<GridX, GridY>;

using DVectorFieldMemXY = VectorFieldMem<double, IdxRangeXY, Direction>;
using DVectorConstFieldXY = VectorField<double const, IdxRangeXY, Direction>;
using DVectorConstFieldSliceXY
        = VectorField<double const, IdxRangeXY, Direction, std::experimental::layout_stride>;

using IdxYX = Idx<GridY, GridX>;
using IdxStepYX = IdxStep<GridY, GridX>;
using IdxRangeYX = IdxRange<GridY, GridX>;

using DVectorFieldMemYX = VectorFieldMem<double, IdxRangeYX, Direction>;
using DVectorConstFieldYX = VectorField<double const, IdxRangeYX, Direction>;


static Idx0D constexpr lbound_0d {};
static IdxStep0D constexpr nelems_0d {};
static IdxRange0D constexpr idx_range_0d(lbound_0d, nelems_0d);

static IdxX constexpr lbound_x(50);
static IdxStepX constexpr nelems_x(3);
static IdxRangeX constexpr idx_range_x(lbound_x, nelems_x);

static IdxY constexpr lbound_y(4);
static IdxStepY constexpr nelems_y(12);

static IdxXY constexpr lbound_x_y {lbound_x, lbound_y};
static IdxStepXY constexpr nelems_x_y(nelems_x, nelems_y);
static IdxRangeXY constexpr idx_range_x_y(lbound_x_y, nelems_x_y);

} // namespace

// Member types of VectorFieldMem 1D \{

TEST(VectorFieldMem0DTest, LayoutType)
{
    EXPECT_TRUE((std::is_same_v<
                 DVectorFieldMem0D::chunk_type::layout_type,
                 std::experimental::layout_right>));
}

TEST(VectorField1DTest, LayoutType)
{
    DVectorFieldMemX field(idx_range_x);

    EXPECT_TRUE((std::is_same_v<
                 std::decay_t<decltype(field)::chunk_type>::layout_type,
                 std::experimental::layout_right>));
}

// TODO: many missing types

// \}
// Functions implemented in VectorFieldMem 1D (and free functions specific to it) \{

TEST(VectorFieldMem0DTest, MoveConstructor)
{
    Coord2D constexpr factor(1.391, 2.444);
    DVectorFieldMem0D field(idx_range_0d);
    ddcHelper::get<Tag1>(field)() = ddc::get<Tag1>(factor);
    ddcHelper::get<Tag2>(field)() = ddc::get<Tag2>(factor);

    DVectorFieldMem0D field2(std::move(field));
    EXPECT_EQ(get_idx_range(field2), idx_range_0d);
    EXPECT_DOUBLE_EQ(ddc::get<Tag1>(factor), ddc::get<Tag1>(field2()));
    EXPECT_DOUBLE_EQ(ddc::get<Tag2>(factor), ddc::get<Tag2>(field2()));
}

TEST(VectorField1DTest, MoveConstructor)
{
    Coord2D constexpr factor(1.391, 2.444);
    DVectorFieldMemX field(idx_range_x);
    for (IdxX ix : idx_range_x) {
        Coord2D val = double((ix - idx_range_x.front()).value()) * factor;
        ddcHelper::get<Tag1>(field)(ix) = ddc::get<Tag1>(val);
        ddcHelper::get<Tag2>(field)(ix) = ddc::get<Tag2>(val);
    }

    DVectorFieldMemX field2(std::move(field));
    EXPECT_EQ(get_idx_range(field2), idx_range_x);
    for (IdxX ix : get_idx_range(field2)) {
        // we expect exact equality, not EXPECT_DOUBLE_EQ: this is the same ref twice
        Coord2D val = double((ix - idx_range_x.front()).value()) * factor;
        EXPECT_EQ(ddc::get<Tag1>(val), ddc::get<Tag1>(field2(ix)));
        EXPECT_EQ(ddc::get<Tag2>(val), ddc::get<Tag2>(field2(ix)));
    }
}

TEST(VectorFieldMem0DTest, MoveAssignment)
{
    Coord2D constexpr factor(1.391, 2.444);
    DVectorFieldMem0D field(idx_range_0d);
    ddcHelper::get<Tag1>(field)() = ddc::get<Tag1>(factor);
    ddcHelper::get<Tag2>(field)() = ddc::get<Tag2>(factor);

    DVectorFieldMem0D field2(IdxRange0D(lbound_0d, IdxStep0D()));
    field2 = std::move(field);
    EXPECT_EQ(get_idx_range(field2), idx_range_0d);
    EXPECT_DOUBLE_EQ(ddc::get<Tag1>(factor), ddc::get<Tag1>(field2()));
    EXPECT_DOUBLE_EQ(ddc::get<Tag2>(factor), ddc::get<Tag2>(field2()));
}

TEST(VectorField1DTest, MoveAssignment)
{
    Coord2D constexpr factor(1.391, 2.444);
    DVectorFieldMemX field(idx_range_x);
    for (IdxX ix : get_idx_range(field)) {
        Coord2D val = double((ix - idx_range_x.front()).value()) * factor;
        ddcHelper::get<Tag1>(field)(ix) = ddc::get<Tag1>(val);
        ddcHelper::get<Tag2>(field)(ix) = ddc::get<Tag2>(val);
    }

    DVectorFieldMemX field2(IdxRangeX(lbound_x, IdxStepX(0)));
    field2 = std::move(field);
    EXPECT_EQ(get_idx_range(field2), idx_range_x);
    for (IdxX ix : get_idx_range(field2)) {
        // we expect exact equality, not EXPECT_DOUBLE_EQ: this is the same ref twice
        Coord2D val = double((ix - idx_range_x.front()).value()) * factor;
        EXPECT_EQ(ddc::get<Tag1>(val), ddc::get<Tag1>(field2(ix)));
        EXPECT_EQ(ddc::get<Tag2>(val), ddc::get<Tag2>(field2(ix)));
    }
}

TEST(VectorFieldMem0DTest, Swap)
{
    Coord2D constexpr factor(1.391, 2.444);
    DVectorFieldMem0D field(idx_range_0d);
    ddcHelper::get<Tag1>(field)() = ddc::get<Tag1>(factor);
    ddcHelper::get<Tag2>(field)() = ddc::get<Tag2>(factor);

    IdxRange0D empty_idx_range(lbound_0d, IdxStep0D());
    DVectorFieldMem0D field2(empty_idx_range);

    std::swap(field2, field);
    EXPECT_EQ(get_idx_range(field), empty_idx_range);
    EXPECT_EQ(get_idx_range(field2), idx_range_0d);
    EXPECT_DOUBLE_EQ(ddc::get<Tag1>(factor), ddc::get<Tag1>(field2()));
    EXPECT_DOUBLE_EQ(ddc::get<Tag2>(factor), ddc::get<Tag2>(field2()));
}

TEST(VectorField1DTest, Swap)
{
    Coord2D constexpr factor(1.391, 2.444);
    DVectorFieldMemX field(idx_range_x);
    for (IdxX ix : get_idx_range(field)) {
        Coord2D val = double((ix - idx_range_x.front()).value()) * factor;
        ddcHelper::get<Tag1>(field)(ix) = ddc::get<Tag1>(val);
        ddcHelper::get<Tag2>(field)(ix) = ddc::get<Tag2>(val);
    }

    IdxRangeX empty_idx_range(lbound_x, IdxStepX(0));
    DVectorFieldMemX field2(empty_idx_range);

    std::swap(field2, field);
    EXPECT_EQ(get_idx_range(field), empty_idx_range);
    EXPECT_EQ(get_idx_range(field2), idx_range_x);
    for (IdxX ix : get_idx_range(field2)) {
        // we expect exact equality, not EXPECT_DOUBLE_EQ: this is the same ref twice
        Coord2D val = double((ix - idx_range_x.front()).value()) * factor;
        EXPECT_EQ(ddc::get<Tag1>(val), ddc::get<Tag1>(field2(ix)));
        EXPECT_EQ(ddc::get<Tag2>(val), ddc::get<Tag2>(field2(ix)));
    }
}

// no dim subset access in 1D

TEST(VectorField1DTest, AccessConst)
{
    Coord2D constexpr factor(1.391, 2.444);
    DVectorFieldMemX field(idx_range_x);
    DVectorFieldMemX const& field_cref = field;
    for (IdxX ix : get_idx_range(field)) {
        Coord2D val = double((ix - idx_range_x.front()).value()) * factor;
        ddcHelper::get<Tag1>(field)(ix) = ddc::get<Tag1>(val);
        ddcHelper::get<Tag2>(field)(ix) = ddc::get<Tag2>(val);
        // we expect exact equality, not EXPECT_DOUBLE_EQ: this is the same ref twice
        EXPECT_EQ(ddc::get<Tag1>(val), ddc::get<Tag1>(field_cref(ix)));
        EXPECT_EQ(ddc::get<Tag2>(val), ddc::get<Tag2>(field_cref(ix)));
    }
}

TEST(VectorField1DTest, Access)
{
    Coord2D constexpr factor(1.391, 2.444);
    DVectorFieldMemX field(idx_range_x);
    for (IdxX ix : get_idx_range(field)) {
        Coord2D val = double((ix - idx_range_x.front()).value()) * factor;
        ddcHelper::get<Tag1>(field)(ix) = ddc::get<Tag1>(val);
        ddcHelper::get<Tag2>(field)(ix) = ddc::get<Tag2>(val);
        // we expect exact equality, not EXPECT_DOUBLE_EQ: this is the same ref twice
        EXPECT_EQ(ddc::get<Tag1>(val), ddc::get<Tag1>(field(ix)));
        EXPECT_EQ(ddc::get<Tag2>(val), ddc::get<Tag2>(field(ix)));
    }
}

TEST(VectorField1DTest, GetConstField)
{
    Coord2D constexpr factor(1.391, 2.444);
    DVectorFieldMemX field(idx_range_x);
    for (IdxX ix : get_idx_range(field)) {
        Coord2D val = double((ix - idx_range_x.front()).value()) * factor;
        ddcHelper::get<Tag1>(field)(ix) = ddc::get<Tag1>(val);
        ddcHelper::get<Tag2>(field)(ix) = ddc::get<Tag2>(val);
        // we expect exact equality, not EXPECT_DOUBLE_EQ: this is the same ref twice
        EXPECT_EQ(ddc::get<Tag1>(val), ddc::get<Tag1>(get_const_field(field)(ix)));
        EXPECT_EQ(ddc::get<Tag2>(val), ddc::get<Tag2>(get_const_field(field)(ix)));
    }
}

TEST(VectorField1DTest, GetFieldFromConst)
{
    Coord2D constexpr factor(1.391, 2.444);
    DVectorFieldMemX field(idx_range_x);
    DVectorFieldMemX const& field_cref = field;
    for (IdxX ix : get_idx_range(field)) {
        Coord2D val = double((ix - idx_range_x.front()).value()) * factor;
        ddcHelper::get<Tag1>(field)(ix) = ddc::get<Tag1>(val);
        ddcHelper::get<Tag2>(field)(ix) = ddc::get<Tag2>(val);
        // we expect exact equality, not EXPECT_DOUBLE_EQ: this is the same ref twice
        EXPECT_EQ(ddc::get<Tag1>(val), ddc::get<Tag1>(get_field(field_cref)(ix)));
        EXPECT_EQ(ddc::get<Tag2>(val), ddc::get<Tag2>(get_field(field_cref)(ix)));
    }
}

TEST(VectorField1DTest, GetField)
{
    Coord2D constexpr factor(1.391, 2.444);
    DVectorFieldMemX field_alloc(idx_range_x);
    DVectorFieldX field = get_field(field_alloc);
    for (IdxX ix : get_idx_range(field)) {
        Coord2D val = double((ix - idx_range_x.front()).value()) * factor;
        ddcHelper::get<Tag1>(field)(ix) = ddc::get<Tag1>(val);
        ddcHelper::get<Tag2>(field)(ix) = ddc::get<Tag2>(val);
        // we expect exact equality, not EXPECT_DOUBLE_EQ: this is the same ref twice
        EXPECT_EQ(ddc::get<Tag1>(val), ddc::get<Tag1>(field_alloc(ix)));
        EXPECT_EQ(ddc::get<Tag2>(val), ddc::get<Tag2>(field_alloc(ix)));
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
    DVectorFieldMemX field(idx_range_x);
    EXPECT_EQ(DVectorFieldMemX::rank(), 1);
}

// TODO: stride

// swap is hidden

TEST(VectorField1DTest, IdxRange)
{
    DVectorFieldMemX field(idx_range_x);
    EXPECT_EQ(idx_range_x, get_idx_range(field));
}

TEST(VectorField1DTest, IdxRangeX)
{
    DVectorFieldMemX field(idx_range_x);
    EXPECT_EQ(idx_range_x, get_idx_range<GridX>(field));
}

// TODO: data_handle()

// TODO: internal_mdspan

// TODO: allocation_mdspan


TEST(VectorField1DTest, Deepcopy)
{
    DVectorFieldMemX field(idx_range_x);
    for (IdxX ix : get_idx_range(field)) {
        ddcHelper::get<Tag1>(field)(ix) = 1.001 * (ix - idx_range_x.front()).value();
        ddcHelper::get<Tag2>(field)(ix) = 0.0;
    }
    DVectorFieldMemX field2(get_idx_range(field));
    ddcHelper::deepcopy(field2, field);
    for (IdxX ix : get_idx_range(field)) {
        // we expect exact equality, not EXPECT_DOUBLE_EQ: these are copy
        EXPECT_EQ(ddc::get<Tag1>(field2(ix)), ddc::get<Tag1>(field(ix)));
        EXPECT_EQ(ddc::get<Tag2>(field2(ix)), ddc::get<Tag2>(field(ix)));
    }
}

// \}
// Functions implemented in VectorFieldMem 2D (and free functions specific to it) \{

// TODO: lots to do still!

TEST(VectorField2DTest, Access)
{
    DVectorFieldMemXY field(idx_range_x_y);
    for (IdxX ix : get_idx_range<GridX>(field)) {
        for (IdxY iy : get_idx_range<GridY>(field)) {
            ddcHelper::get<Tag1>(field)(ix, iy) = 1.357 * (ix - idx_range_x.front()).value();
            ddcHelper::get<Tag2>(field)(ix, iy)
                    = 1.159 * (iy - ddc::select<GridY>(idx_range_x_y).front()).value();
            // we expect exact equality, not EXPECT_DOUBLE_EQ: this is the same ref twice
            EXPECT_EQ(ddc::get<Tag1>(field(ix, iy)), ddc::get<Tag1>(field(ix, iy)));
            EXPECT_EQ(ddc::get<Tag2>(field(ix, iy)), ddc::get<Tag2>(field(ix, iy)));
        }
    }
}

TEST(VectorField2DTest, AccessReordered)
{
    DVectorFieldMemXY field(idx_range_x_y);
    for (IdxX ix : get_idx_range<GridX>(field)) {
        for (IdxY iy : get_idx_range<GridY>(field)) {
            ddcHelper::get<Tag1>(field)(ix, iy) = 1.455 * (ix - idx_range_x.front()).value();
            ddcHelper::get<Tag2>(field)(ix, iy)
                    = 1.522 * (iy - ddc::select<GridY>(idx_range_x_y).front()).value();
            // we expect exact equality, not EXPECT_DOUBLE_EQ: this is the same ref twice
            EXPECT_EQ(ddc::get<Tag1>(field(iy, ix)), ddc::get<Tag1>(field(ix, iy)));
            EXPECT_EQ(ddc::get<Tag2>(field(iy, ix)), ddc::get<Tag2>(field(ix, iy)));
        }
    }
}

TEST(VectorField2DTest, Cview)
{
    DVectorFieldMemXY field(idx_range_x_y);
    for (IdxX ix : get_idx_range<GridX>(field)) {
        for (IdxY iy : get_idx_range<GridY>(field)) {
            ddcHelper::get<Tag1>(field)(ix, iy) = 1. * (ix - idx_range_x.front()).value();
            ddcHelper::get<Tag2>(field)(ix, iy)
                    = .001 * (iy - ddc::select<GridY>(idx_range_x_y).front()).value();
        }
    }
    DVectorConstFieldXY cview = get_const_field(field);
    for (IdxX ix : get_idx_range<GridX>(field)) {
        for (IdxY iy : get_idx_range<GridY>(field)) {
            // we expect complete equality, not EXPECT_DOUBLE_EQ: these are copy
            EXPECT_EQ(ddc::get<Tag1>(cview(ix, iy)), ddc::get<Tag1>(field(ix, iy)));
            EXPECT_EQ(ddc::get<Tag2>(cview(ix, iy)), ddc::get<Tag2>(field(ix, iy)));
        }
    }
}

TEST(VectorField2DTest, SliceCoordX)
{
    IdxX constexpr slice_x_val = IdxX(lbound_x + 1);

    DVectorFieldMemXY field(idx_range_x_y);
    DVectorFieldMemXY const& field_cref = field;
    for (IdxX ix : get_idx_range<GridX>(field)) {
        for (IdxY iy : get_idx_range<GridY>(field)) {
            ddcHelper::get<Tag1>(field)(ix, iy) = 1. * (ix - idx_range_x.front()).value();
            ddcHelper::get<Tag2>(field)(ix, iy)
                    = .001 * (iy - ddc::select<GridY>(idx_range_x_y).front()).value();
        }
    }

    DVectorConstFieldY field_y = field_cref[slice_x_val];
    EXPECT_TRUE((std::is_same_v<
                 std::decay_t<decltype(field_y)>::layout_type,
                 std::experimental::layout_right>));
    EXPECT_EQ(
            ddcHelper::get<Tag1>(field_y).extent<GridY>(),
            ddcHelper::get<Tag1>(field).extent<GridY>());
    EXPECT_EQ(
            ddcHelper::get<Tag2>(field_y).extent<GridY>(),
            ddcHelper::get<Tag2>(field).extent<GridY>());
    for (IdxY ix : get_idx_range<GridY>(field_cref)) {
        // we expect complete equality, not EXPECT_DOUBLE_EQ: these are copy
        EXPECT_EQ(ddc::get<Tag1>(field_y(ix)), ddc::get<Tag1>(field_cref(slice_x_val, ix)));
        EXPECT_EQ(ddc::get<Tag2>(field_y(ix)), ddc::get<Tag2>(field_cref(slice_x_val, ix)));
    }
}

TEST(VectorField2DTest, SliceCoordY)
{
    IdxY constexpr slice_y_val = IdxY(lbound_y + 1);

    DVectorFieldMemXY field(idx_range_x_y);
    DVectorFieldMemXY const& field_cref = field;
    for (IdxX ix : get_idx_range<GridX>(field)) {
        for (IdxY iy : get_idx_range<GridY>(field)) {
            ddcHelper::get<Tag1>(field)(ix, iy) = 1. * (ix - idx_range_x.front()).value();
            ddcHelper::get<Tag2>(field)(ix, iy)
                    = .001 * (iy - ddc::select<GridY>(idx_range_x_y).front()).value();
        }
    }

    DVectorConstFieldSliceX field_x = field_cref[slice_y_val];
    EXPECT_TRUE((std::is_same_v<
                 std::decay_t<decltype(field_x)>::layout_type,
                 std::experimental::layout_stride>));
    EXPECT_EQ(
            ddcHelper::get<Tag1>(field_x).extent<GridX>(),
            ddcHelper::get<Tag1>(field).extent<GridX>());
    EXPECT_EQ(
            ddcHelper::get<Tag2>(field_x).extent<GridX>(),
            ddcHelper::get<Tag2>(field).extent<GridX>());
    for (IdxX ix : get_idx_range<GridX>(field_cref)) {
        // we expect complete equality, not EXPECT_DOUBLE_EQ: these are copy
        EXPECT_EQ(ddc::get<Tag1>(field_x(ix)), ddc::get<Tag1>(field_cref(ix, slice_y_val)));
        EXPECT_EQ(ddc::get<Tag2>(field_x(ix)), ddc::get<Tag2>(field_cref(ix, slice_y_val)));
    }
}

TEST(VectorField2DTest, IdxRangeSliceX)
{
    IdxRangeX constexpr subidx_range_x = IdxRangeX(IdxX(lbound_x + 1), IdxStepX(nelems_x - 2));

    DVectorFieldMemXY field(idx_range_x_y);
    DVectorFieldMemXY const& field_cref = field;
    for (IdxX ix : get_idx_range<GridX>(field)) {
        for (IdxY iy : get_idx_range<GridY>(field)) {
            ddcHelper::get<Tag1>(field)(ix, iy) = 1. * (ix - idx_range_x.front()).value();
            ddcHelper::get<Tag2>(field)(ix, iy)
                    = .001 * (iy - ddc::select<GridY>(idx_range_x_y).front()).value();
        }
    }

    DVectorConstFieldXY subfield_x = field_cref[subidx_range_x];
    EXPECT_TRUE((std::is_same_v<
                 std::decay_t<decltype(subfield_x)>::layout_type,
                 std::experimental::layout_right>));

    EXPECT_EQ(ddcHelper::get<Tag1>(subfield_x).extent<GridX>(), subidx_range_x.size());
    EXPECT_EQ(ddcHelper::get<Tag2>(subfield_x).extent<GridY>(), get_idx_range<GridY>(field).size());
    EXPECT_EQ(ddcHelper::get<Tag1>(subfield_x).extent<GridX>(), subidx_range_x.size());
    EXPECT_EQ(ddcHelper::get<Tag2>(subfield_x).extent<GridY>(), get_idx_range<GridY>(field).size());
    for (IdxX ix : get_idx_range<GridX>(subfield_x)) {
        for (IdxY iy : get_idx_range<GridY>(subfield_x)) {
            // we expect complete equality, not EXPECT_DOUBLE_EQ: these are copy
            EXPECT_EQ(ddc::get<Tag1>(subfield_x(ix, iy)), ddc::get<Tag1>(field_cref(ix, iy)));
            EXPECT_EQ(ddc::get<Tag2>(subfield_x(ix, iy)), ddc::get<Tag2>(field_cref(ix, iy)));
        }
    }
}

TEST(VectorField2DTest, IdxRangeSliceXTooearly)
{
    [[maybe_unused]] IdxRangeX constexpr subidx_range_x = IdxRangeX(IdxX(lbound_x - 1), nelems_x);

    DVectorFieldMemXY field(idx_range_x_y);
#ifndef NDEBUG // The assertion is only checked if NDEBUG isn't defined
    // the error message is checked with clang & gcc only
    EXPECT_DEATH(
            field[subidx_range_x],
            R"rgx([Aa]ssert.*uid<ODDims>\(m_element_begin\).*uid<ODDims>\(odomain\.m_element_begin\))rgx");
#endif
}

TEST(VectorField2DTest, IdxRangeSliceXToolate)
{
    [[maybe_unused]] IdxRangeX constexpr subidx_range_x
            = IdxRangeX(lbound_x, IdxStepX(nelems_x + 1));

    DVectorFieldMemXY field(idx_range_x_y);
#ifndef NDEBUG // The assertion is only checked if NDEBUG isn't defined
    // the error message is checked with clang & gcc only
    EXPECT_DEATH(
            field[subidx_range_x],
            R"rgx([Aa]ssert.*uid<ODDims>\(m_element_end\).*uid<ODDims>\(odomain\.m_element_end\).*)rgx");
#endif
}

TEST(VectorField2DTest, IdxRangeSliceY)
{
    IdxRangeY constexpr subidx_range_y = IdxRangeY(IdxY(lbound_y + 1), IdxStepY(nelems_y - 2));

    DVectorFieldMemXY field(idx_range_x_y);
    DVectorFieldMemXY const& field_cref = field;
    for (IdxX ix : get_idx_range<GridX>(field)) {
        for (IdxY iy : get_idx_range<GridY>(field)) {
            ddcHelper::get<Tag1>(field)(ix, iy) = 1. * (ix - idx_range_x.front()).value();
            ddcHelper::get<Tag2>(field)(ix, iy)
                    = .001 * (iy - ddc::select<GridY>(idx_range_x_y).front()).value();
        }
    }
    DVectorConstFieldSliceXY subfield_y = field_cref[subidx_range_y];
    EXPECT_TRUE((std::is_same_v<
                 std::decay_t<decltype(subfield_y)>::layout_type,
                 std::experimental::layout_stride>));

    EXPECT_EQ(ddcHelper::get<Tag1>(subfield_y).extent<GridX>(), get_idx_range<GridX>(field).size());
    EXPECT_EQ(ddcHelper::get<Tag1>(subfield_y).extent<GridY>(), subidx_range_y.size());
    EXPECT_EQ(ddcHelper::get<Tag2>(subfield_y).extent<GridX>(), get_idx_range<GridX>(field).size());
    EXPECT_EQ(ddcHelper::get<Tag2>(subfield_y).extent<GridY>(), subidx_range_y.size());
    for (IdxX ix : get_idx_range<GridX>(subfield_y)) {
        for (IdxY iy : get_idx_range<GridY>(subfield_y)) {
            // we expect complete equality, not EXPECT_DOUBLE_EQ: these are copy
            EXPECT_EQ(ddc::get<Tag1>(subfield_y(ix, iy)), ddc::get<Tag1>(field_cref(ix, iy)));
            EXPECT_EQ(ddc::get<Tag2>(subfield_y(ix, iy)), ddc::get<Tag2>(field_cref(ix, iy)));
        }
    }
}

TEST(VectorField2DTest, Deepcopy)
{
    DVectorFieldMemXY field(idx_range_x_y);
    for (IdxX ix : get_idx_range<GridX>(field)) {
        for (IdxY iy : get_idx_range<GridY>(field)) {
            ddcHelper::get<Tag1>(field)(ix, iy) = 1.739 * (ix - idx_range_x.front()).value();
            ddcHelper::get<Tag2>(field)(ix, iy)
                    = 1.412 * (iy - ddc::select<GridY>(idx_range_x_y).front()).value();
        }
    }
    DVectorFieldMemXY field2(get_idx_range(field));
    ddcHelper::deepcopy(field2, field);
    for (IdxX ix : get_idx_range<GridX>(field)) {
        for (IdxY iy : get_idx_range<GridY>(field)) {
            // we expect complete equality, not EXPECT_DOUBLE_EQ: these are copy
            EXPECT_EQ(ddc::get<Tag1>(field2(ix, iy)), ddc::get<Tag1>(field(ix, iy)));
            EXPECT_EQ(ddc::get<Tag2>(field2(ix, iy)), ddc::get<Tag2>(field(ix, iy)));
        }
    }
}

TEST(VectorField2DTest, DeepcopyReordered)
{
    DVectorFieldMemXY field(idx_range_x_y);
    for (IdxX ix : get_idx_range<GridX>(field)) {
        for (IdxY iy : get_idx_range<GridY>(field)) {
            ddcHelper::get<Tag1>(field)(ix, iy) = 1.739 * (ix - idx_range_x.front()).value();
            ddcHelper::get<Tag2>(field)(ix, iy)
                    = 1.412 * (iy - ddc::select<GridY>(idx_range_x_y).front()).value();
        }
    }
    DVectorFieldMemYX field2_alloc(ddc::select<GridY, GridX>(get_idx_range(field)));
    VectorField<double, IdxRangeXY, Direction, std::experimental::layout_left>
            field2(get_idx_range(field),
                   field2_alloc.get<Tag1>().data_handle(),
                   field2_alloc.get<Tag2>().data_handle());
    ddcHelper::deepcopy(field2, field);
    for (IdxX ix : get_idx_range<GridX>(field)) {
        for (IdxY iy : get_idx_range<GridY>(field)) {
            // we expect complete equality, not EXPECT_DOUBLE_EQ: these are copy
            EXPECT_EQ(ddc::get<Tag1>(field2_alloc(ix, iy)), ddc::get<Tag1>(field(ix, iy)));
            EXPECT_EQ(ddc::get<Tag2>(field2_alloc(ix, iy)), ddc::get<Tag2>(field(ix, iy)));
            EXPECT_EQ(ddc::get<Tag1>(field2_alloc(ix, iy)), ddc::get<Tag1>(field(iy, ix)));
            EXPECT_EQ(ddc::get<Tag2>(field2_alloc(ix, iy)), ddc::get<Tag2>(field(iy, ix)));
        }
    }
}
