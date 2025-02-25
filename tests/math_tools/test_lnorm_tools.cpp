// SPDX-License-Identifier: MIT
#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "directional_tag.hpp"
#include "l_norm_tools.hpp"
#include "vector_field.hpp"
#include "vector_field_mem.hpp"

namespace {

struct X
{
    /// A boolean indicating if dimension describes a covariant coordinate.
    static bool constexpr IS_COVARIANT = true;
    /// A boolean indicating if dimension describes a contravariant coordinate.
    static bool constexpr IS_CONTRAVARIANT = true;
    /// A type-alias mapping to the co/contra-variant counterpart.
    using Dual = X;
};
struct Y
{
    /// A boolean indicating if dimension describes a covariant coordinate.
    static bool constexpr IS_COVARIANT = true;
    /// A boolean indicating if dimension describes a contravariant coordinate.
    static bool constexpr IS_CONTRAVARIANT = true;
    /// A type-alias mapping to the co/contra-variant counterpart.
    using Dual = Y;
};

using GridX = UniformGridBase<X>;
using GridY = UniformGridBase<Y>;

using IdxX = Idx<GridX>;
using IdxY = Idx<GridY>;
using IdxXY = Idx<GridX, GridY>;

using IdxStepX = IdxStep<GridX>;
using IdxStepY = IdxStep<GridY>;
using IdxStepXY = IdxStep<GridX, GridY>;

using IdxRangeX = IdxRange<GridX>;
using IdxRangeY = IdxRange<GridY>;
using IdxRangeXY = IdxRange<GridX, GridY>;

using CoordX = Coord<X>;
using CoordY = Coord<Y>;
using CoordXY = Coord<X, Y>;

using DFieldMemXY = DFieldMem<IdxRangeXY>;
using DVectorFieldMemXY = VectorFieldMem<double, IdxRangeXY, NDTag<X, Y>>;

using DFieldXY = DField<IdxRangeXY>;
using DVectorFieldXY = VectorField<double, IdxRangeXY, NDTag<X, Y>>;

template <class Grid1D>
KOKKOS_FUNCTION typename Grid1D::continuous_element_type get_coordinate(Idx<Grid1D> x)
{
    using Dim = typename Grid1D::continuous_dimension_type;
    CoordXY const origin(0.0, 0.0);
    CoordXY step(0.2, 0.1);
    Idx<Grid1D> const origin_idx(0);

    return ddc::select<Dim>(origin) + ddc::select<Dim>(step) * (x - origin_idx);
}

void test_norm_inf_field()
{
    IdxX idx_start_x(0);
    IdxY idx_start_y(2);
    IdxStepX idx_step_x(5);
    IdxStepY idx_step_y(10);
    IdxRangeX idx_range_x(idx_start_x, idx_step_x);
    IdxRangeY idx_range_y(idx_start_y, idx_step_y);
    IdxRangeXY idx_range_xy(idx_range_x, idx_range_y);

    DFieldMemXY field_alloc(idx_range_xy);
    DFieldXY field(field_alloc);
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            idx_range_xy,
            KOKKOS_LAMBDA(IdxXY idx) {
                const double x = get_coordinate(ddc::select<GridX>(idx));
                const double y = get_coordinate(ddc::select<GridY>(idx));
                field(idx) = 10 - (x - 0.4) * (x - 0.4) * (y - 2.5) * (y - 2.5);
            });
    double max_val = norm_inf(Kokkos::DefaultExecutionSpace(), get_const_field(field));
    EXPECT_NEAR(max_val, 10.0, 1e-14);
}

void test_norm_inf_vector_field()
{
    IdxX idx_start_x(0);
    IdxY idx_start_y(2);
    IdxStepX idx_step_x(5);
    IdxStepY idx_step_y(10);
    IdxRangeX idx_range_x(idx_start_x, idx_step_x);
    IdxRangeY idx_range_y(idx_start_y, idx_step_y);
    IdxRangeXY idx_range_xy(idx_range_x, idx_range_y);

    DVectorFieldMemXY field_alloc(idx_range_xy);
    DVectorFieldXY field(field_alloc);
    DFieldXY field_x(field.get<X>());
    DFieldXY field_y(field.get<Y>());
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            idx_range_xy,
            KOKKOS_LAMBDA(IdxXY idx) {
                const double x = get_coordinate(ddc::select<GridX>(idx));
                const double y = get_coordinate(ddc::select<GridY>(idx));
                field_x(idx) = 10 - (x - 0.4) * (x - 0.4);
                field_y(idx) = 8 - (y - 2.5) * (x - 2.5);
            });
    double max_val = norm_inf(Kokkos::DefaultExecutionSpace(), get_const_field(field));
    EXPECT_NEAR(max_val, 10.0, 1e-14);
}

TEST(MathTools, NormInfField)
{
    test_norm_inf_field();
}

TEST(MathTools, NormInfVectorField)
{
    test_norm_inf_vector_field();
}

} // namespace
