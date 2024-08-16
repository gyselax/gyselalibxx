// SPDX-License-Identifier: MIT
#include <ddc/ddc.hpp>
#include <ddc/kernels/splines/deriv.hpp>

#include <sll/view.hpp>

#include <gtest/gtest.h>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "derivative_field.hpp"
#include "derivative_field_mem.hpp"
#include "directional_tag.hpp"
#include "grid_builder.hpp"

namespace {

struct GridX
{
};
using dX = ddc::Deriv<GridX>;
using IdxX = Idx<GridX>;
using IdxStepX = IdxStep<GridX>;
using IdxRangeX = IdxRange<GridX>;


struct GridY
{
};
using dY = ddc::Deriv<GridY>;
using IdxY = Idx<GridY>;
using IdxStepY = IdxStep<GridY>;
using IdxRangeY = IdxRange<GridY>;

using IdxXY = Idx<GridX, GridY>;
using IdxStepXY = IdxStep<GridX, GridY>;
using IdxRangeXY = IdxRange<GridX, GridY>;

using IdxXdX = Idx<dX, GridX>;

using Idx_dXdYXY = Idx<dX, dY, GridX, GridY>;
using IdxRange_dXdYXY = IdxRange<dX, dY, GridX, GridY>;

static IdxX constexpr lbound_x(50);
static IdxStepX constexpr nelems_x(3);
static IdxRangeX constexpr idx_range_x(lbound_x, nelems_x);

static IdxY constexpr lbound_y(4);
static IdxStepY constexpr nelems_y(12);
static IdxRangeY constexpr idx_range_y(lbound_y, nelems_y);

static Idx<dX> constexpr lbound_dx(0);
static IdxStep<dX> constexpr nelems_dx(4);
static IdxRange<dX> constexpr idx_range_dx(lbound_dx, nelems_dx);

static Idx<dY> constexpr lbound_dy(0);
static IdxStep<dY> constexpr nelems_dy(4);
static IdxRange<dY> constexpr idx_range_dy(lbound_dy, nelems_dy);

static IdxXY constexpr lbound_x_y {lbound_x, lbound_y};
static IdxStepXY constexpr nelems_x_y(nelems_x, nelems_y);
static IdxRangeXY constexpr idx_range_x_y(lbound_x_y, nelems_x_y);

static IdxRange_dXdYXY constexpr idx_range_x_y_dx_dy(idx_range_x_y, idx_range_dx, idx_range_dy);
} // namespace

// Test the constructor for a 2D field with derivatives in 1 direction
TEST(DerivFieldMemTest, Constructor1Deriv)
{
    // Type for a x,y field with 1 derivative in x
    using DFieldMemXY_dX = DerivFieldMem<double, IdxRange<dX, GridX, GridY>, 1>;

    // Index range where derivatives are defined
    IdxRangeSlice<GridX> deriv_idx_range_x(idx_range_x.front(), IdxStepX(2), idx_range_x.extents());

    // Define the field memory allocation
    DFieldMemXY_dX dxField(idx_range_x_y, deriv_idx_range_x);

    // Ensure that the internal field has the expected type
    constexpr bool same = std::is_same_v<
            typename DFieldMemXY_dX::chunk_type,
            host_t<DFieldMem<IdxRange<dX, GridX, GridY>>>>;
    EXPECT_TRUE(same);
}

// Test the constructor for a 2D field with derivatives in 2 directions
TEST(DerivFieldMemTest, Constructor2Deriv)
{
    // Type for a x,y field with 1 derivative in x and 1 derivative in y
    using DFieldMemXY_dXdY = DerivFieldMem<double, IdxRange_dXdYXY, 1>;

    // Index ranges where derivatives are defined
    IdxRangeSlice<GridX> deriv_idx_range_x(idx_range_x.front(), IdxStepX(2), idx_range_x.extents());
    IdxRangeSlice<GridY> deriv_idx_range_y(idx_range_y.front(), IdxStepY(2), idx_range_y.extents());

    // Define the field memory allocation
    DFieldMemXY_dXdY dxdyField(idx_range_x_y, deriv_idx_range_x, deriv_idx_range_y);

    // Ensure that the internal field has the expected type
    bool same = std::
            is_same_v<typename DFieldMemXY_dXdY::chunk_type, host_t<DFieldMem<IdxRange_dXdYXY>>>;
    EXPECT_TRUE(same);
}

// Test the field constructor for a 2D field with derivatives in 1 direction
TEST(DerivFieldTest, Constructor1Deriv)
{
    // Type for a x,y field with 1 derivative in x
    using DFieldXY_dX = DerivField<double, IdxRange<dX, GridX, GridY>>;

    // Index range where derivatives are defined
    IdxRangeSlice<GridX> deriv_idx_range_x(idx_range_x.front(), IdxStepX(2), idx_range_x.extents());

    // Define a field memory allocation on x-y with 1 derivative in x and y
    DerivFieldMem<double, IdxRange<dX, GridX, GridY>, 1>
            dxField_alloc(idx_range_x_y, deriv_idx_range_x);
    // Define the field
    DFieldXY_dX dxField(dxField_alloc);

    // Ensure that the internal field has the expected type
    bool same = std::is_same_v<DFieldXY_dX::chunk_type, host_t<DField<IdxRange<dX, GridX, GridY>>>>;
    EXPECT_TRUE(same);
}

// Test the field constructor for a 2D field with derivatives in 2 directions
TEST(DerivFieldTest, Constructor2Deriv)
{
    // Type for a x,y field with 1 derivative in x and 1 derivative in y
    using DFieldXY_dXdY = DerivField<double, IdxRange_dXdYXY>;

    // Index ranges where derivatives are defined
    IdxRangeSlice<GridX> deriv_idx_range_x(idx_range_x.front(), IdxStepX(2), idx_range_x.extents());
    IdxRangeSlice<GridY> deriv_idx_range_y(idx_range_y.front(), IdxStepY(2), idx_range_y.extents());

    // Define a field memory allocation on x-y with 1 derivative in x and y
    DerivFieldMem<double, IdxRange_dXdYXY, 1>
            dxdyField_alloc(idx_range_x_y, deriv_idx_range_x, deriv_idx_range_y);
    // Define the field
    DFieldXY_dXdY dxdyField(dxdyField_alloc);

    // Ensure that the internal field has the expected type
    bool same = std::is_same_v<DFieldXY_dXdY::chunk_type, host_t<DField<IdxRange_dXdYXY>>>;
    EXPECT_TRUE(same);
}

// Test if the values of the function can be accessed via the get_values_field function
TEST(DerivFieldTest, FieldValueAccess)
{
    // Index ranges where derivatives are defined
    IdxRangeSlice<GridX> deriv_idx_range_x(idx_range_x.front(), IdxStepX(2), idx_range_x.extents());
    IdxRangeSlice<GridY> deriv_idx_range_y(idx_range_y.front(), IdxStepY(2), idx_range_y.extents());

    // Define a field memory allocationon x-y with 1 derivative in x and y
    DerivFieldMem<double, IdxRange_dXdYXY, 1>
            dxdyField_alloc(idx_range_x_y, deriv_idx_range_x, deriv_idx_range_y);
    // Define the field
    DerivField<double, IdxRange_dXdYXY> dxdyField(dxdyField_alloc);

    // Check that the index range of the values matches the expected index range
    EXPECT_EQ(idx_range_x_y, get_idx_range(dxdyField.get_values_field()));
}

// Test if the values of the function can be accessed via the constant get_values_field function
TEST(DerivFieldTest, ViewValueAccess)
{
    // Index ranges where derivatives are defined
    IdxRangeSlice<GridX> deriv_idx_range_x(idx_range_x.front(), IdxStepX(2), idx_range_x.extents());
    IdxRangeSlice<GridY> deriv_idx_range_y(idx_range_y.front(), IdxStepY(2), idx_range_y.extents());

    // Define a field memory allocation on x-y with 1 derivative in x and y
    DerivFieldMem<double, IdxRange_dXdYXY, 1>
            dxdyField_alloc(idx_range_x_y, deriv_idx_range_x, deriv_idx_range_y);
    // Define the field
    DerivConstField<const double, IdxRange_dXdYXY> dxdyField(dxdyField_alloc);

    // Check that the index range of the values matches the expected index range
    EXPECT_EQ(idx_range_x_y, get_idx_range(dxdyField.get_values_field()));
}


// Test if the derivatives of the function can be accessed via the slice function
TEST(DerivFieldTest, derivValueAccess)
{
    // Index ranges where derivatives are defined
    IdxRangeSlice<GridX> deriv_idx_range_x(idx_range_x.front(), IdxStepX(2), idx_range_x.extents());
    IdxRangeSlice<GridY> deriv_idx_range_y(idx_range_y.front(), IdxStepY(2), idx_range_y.extents());

    // Define a field on x-y with 2 derivatives in x and y
    DerivFieldMem<double, IdxRange_dXdYXY, 2>
            dxdyField(idx_range_x_y, deriv_idx_range_x, deriv_idx_range_y);

    // The derivatives to be retrieved from the field
    IdxRange<dX> x_deriv_block(Idx<dX>(1), IdxStep<dX>(2));
    IdxRange<dY> y_deriv_block(Idx<dY>(1), IdxStep<dY>(2));

    // The subindex ranges to be retrieved from the field
    IdxRange<GridX> deriv_x_subdom(idx_range_x.front(), IdxStep<GridX>(1));
    IdxRange<GridY> deriv_y_subdom(idx_range_y.front(), IdxStep<GridY>(1));

    // The expected index range for the x-derivatives identified in x_deriv_block at the positions identified by deriv_x_subdom
    IdxRange<dX, GridX, GridY>
            dx_idx_range(x_deriv_block, deriv_x_subdom, ddc::select<GridY>(idx_range_x_y));
    // The expected index range for the y-derivatives identified in y_deriv_block at the positions identified by deriv_y_subdom
    IdxRange<dY, GridX, GridY>
            dy_idx_range(y_deriv_block, ddc::select<GridX>(idx_range_x_y), deriv_y_subdom);
    // The expected index range for the x-y-derivatives identified in x_deriv_block and y_deriv_block at the positions identified
    // by deriv_x_subdom and deriv_y_subdom
    IdxRange_dXdYXY dx_dy_idx_range(x_deriv_block, y_deriv_block, deriv_x_subdom, deriv_y_subdom);

    // Build the multi-D index ranges that will slice the field
    IdxRange<dX, GridX> slice_idx_dx(x_deriv_block, deriv_x_subdom);
    IdxRange<dY, GridY> slice_idx_dy(y_deriv_block, deriv_y_subdom);
    IdxRange_dXdYXY slice_idx_dx_dy(x_deriv_block, y_deriv_block, deriv_x_subdom, deriv_y_subdom);

    // Collect the index range of the sliced fields
    IdxRange<dX, GridX, GridY> slice_dx_idx_range = get_idx_range(dxdyField[slice_idx_dx]);
    IdxRange<dY, GridX, GridY> slice_dy_idx_range = get_idx_range(dxdyField[slice_idx_dy]);
    IdxRange_dXdYXY slice_dx_dy_idx_range = get_idx_range(dxdyField[slice_idx_dx_dy]);

    // Check that the index ranges are as expected
    EXPECT_EQ(dx_idx_range, slice_dx_idx_range);
    EXPECT_EQ(dy_idx_range, slice_dy_idx_range);
    EXPECT_EQ(dx_dy_idx_range, slice_dx_dy_idx_range);
}

// Test if the derivatives of the function can be accessed from a field via the slice function
TEST(DerivFieldTest, derivFieldValueAccess)
{
    // Index ranges where derivatives are defined
    IdxRangeSlice<GridX> deriv_idx_range_x(idx_range_x.front(), IdxStepX(2), idx_range_x.extents());
    IdxRangeSlice<GridY> deriv_idx_range_y(idx_range_y.front(), IdxStepY(2), idx_range_y.extents());

    // Define a field memory allocation on x-y with 2 derivatives in x and y
    DerivFieldMem<double, IdxRange_dXdYXY, 2>
            dxdyField_alloc(idx_range_x_y, deriv_idx_range_x, deriv_idx_range_y);
    // Define the field
    DerivField<double, IdxRange_dXdYXY> dxdyField(dxdyField_alloc);

    // The derivatives to be retrieved from the field
    IdxRange<dX> x_deriv_block(Idx<dX>(1), IdxStep<dX>(2));
    IdxRange<dY> y_deriv_block(Idx<dY>(1), IdxStep<dY>(2));

    // The subindex ranges to be retrieved from the field
    IdxRange<GridX> deriv_x_subdom(idx_range_x.front(), IdxStep<GridX>(1));
    IdxRange<GridY> deriv_y_subdom(idx_range_y.front(), IdxStep<GridY>(1));

    // The expected index range for the x-derivatives identified in x_deriv_block at the positions identified by deriv_x_subdom
    IdxRange<dX, GridX, GridY>
            dx_idx_range(x_deriv_block, deriv_x_subdom, ddc::select<GridY>(idx_range_x_y));
    // The expected index range for the y-derivatives identified in y_deriv_block at the positions identified by deriv_y_subdom
    IdxRange<dY, GridX, GridY>
            dy_idx_range(y_deriv_block, ddc::select<GridX>(idx_range_x_y), deriv_y_subdom);
    // The expected index range for the x-y-derivatives identified in x_deriv_block and y_deriv_block at the positions identified
    // by deriv_x_subdom and deriv_y_subdom
    IdxRange_dXdYXY dx_dy_idx_range(x_deriv_block, y_deriv_block, deriv_x_subdom, deriv_y_subdom);

    // Build the multi-D index ranges that will slice the field
    IdxRange<dX, GridX> slice_idx_dx(x_deriv_block, deriv_x_subdom);
    IdxRange<dY, GridY> slice_idx_dy(y_deriv_block, deriv_y_subdom);
    IdxRange_dXdYXY slice_idx_dx_dy(x_deriv_block, y_deriv_block, deriv_x_subdom, deriv_y_subdom);

    // Collect the index range of the sliced fields
    IdxRange<dX, GridX, GridY> slice_dx_idx_range = get_idx_range(dxdyField[slice_idx_dx]);
    IdxRange<dY, GridX, GridY> slice_dy_idx_range = get_idx_range(dxdyField[slice_idx_dy]);
    IdxRange_dXdYXY slice_dx_dy_idx_range = get_idx_range(dxdyField[slice_idx_dx_dy]);

    // Check that the index ranges are as expected
    EXPECT_EQ(dx_idx_range, slice_dx_idx_range);
    EXPECT_EQ(dy_idx_range, slice_dy_idx_range);
    EXPECT_EQ(dx_dy_idx_range, slice_dx_dy_idx_range);
}

// Test the element-wise operator
TEST(DerivFieldTest, ElementAccess)
{
    // Index ranges where derivatives are defined
    IdxRangeSlice<GridX> deriv_idx_range_x(idx_range_x.front(), IdxStepX(2), idx_range_x.extents());
    IdxRangeSlice<GridY> deriv_idx_range_y(idx_range_y.front(), IdxStepY(2), idx_range_y.extents());

    // Define a field memory allocation on x-y with 3 derivatives in x and y
    DerivFieldMem<int, IdxRange_dXdYXY, 3>
            dxdyField_alloc(idx_range_x_y, deriv_idx_range_x, deriv_idx_range_y);
    // Define the field
    DerivField<int, IdxRange_dXdYXY> dxdyField(dxdyField_alloc);

    // A subset  of the x-derivatives to be retrieved with get_mdspan
    IdxRange<dX> x_deriv_block(Idx<dX>(1), IdxStep<dX>(2));
    // A subset  of the y-derivatives to be retrieved with get_mdspan
    IdxRange<dY> y_deriv_block(Idx<dY>(1), IdxStep<dY>(2));
    // A subset  of the x and y derivatives to be retrieved with get_mdspan
    IdxRange<dX, dY> x_y_deriv_block(x_deriv_block, y_deriv_block);

    // Get an mdspan describing the values of the function
    detail::ViewNDMaker<2, int, false>::type vals = dxdyField.get_mdspan();
    // Check the shape of the mdspan
    EXPECT_EQ(vals.rank(), 2);
    EXPECT_EQ(vals.extent(0), idx_range_x.size());
    EXPECT_EQ(vals.extent(1), idx_range_y.size());

    GridBuilder example_grid(idx_range_x_y_dx_dy);

    Idx<dX, dY> no_deriv_idx(0, 0);
    // Fill the mdspan with values
    ddc::for_each(idx_range_x_y, [&](Idx<GridX, GridY> idx_x_y) {
        vals(ddc::select<GridX>(idx_x_y) - idx_range_x.front(),
             ddc::select<GridY>(idx_x_y) - idx_range_y.front())
                = example_grid(Idx_dXdYXY(idx_x_y, no_deriv_idx));
    });

    // Get an mdspan describing the 1st and 2nd x-derivatives of the function
    detail::ViewNDMaker<3, int, false>::type dx_vals = dxdyField.get_mdspan(x_deriv_block);
    // Check the shape of the mdspan
    EXPECT_EQ(dx_vals.rank(), 3);
    EXPECT_EQ(dx_vals.extent(0), x_deriv_block.size());
    EXPECT_EQ(dx_vals.extent(1), deriv_idx_range_x.size());
    EXPECT_EQ(dx_vals.extent(2), idx_range_y.size());

    Idx<dY> no_deriv_idx_y(0);
    // Fill the mdspan with values
    for (IdxX idx_x : deriv_idx_range_x) {
        ddc::for_each(idx_range_y, [&](IdxY idx_y) {
            ddc::for_each(x_deriv_block, [&](Idx<dX> idx_dx) {
                dx_vals(idx_dx - x_deriv_block.front(),
                        deriv_idx_range_x.get_index(idx_x),
                        idx_y - idx_range_y.front())
                        = example_grid(Idx_dXdYXY(idx_x, idx_y, idx_dx, no_deriv_idx_y));
            });
        });
    }

    // Get an mdspan describing the 1st and 2nd y-derivatives of the function
    detail::ViewNDMaker<3, int, false>::type dy_vals = dxdyField.get_mdspan(y_deriv_block);
    // Check the shape of the mdspan
    EXPECT_EQ(dy_vals.rank(), 3);
    EXPECT_EQ(dy_vals.extent(0), y_deriv_block.size());
    EXPECT_EQ(dy_vals.extent(1), idx_range_x.size());
    EXPECT_EQ(dy_vals.extent(2), deriv_idx_range_y.size());

    Idx<dX> no_deriv_idx_x(0);
    // Fill the mdspan with values
    for (IdxY idx_y : deriv_idx_range_y) {
        ddc::for_each(idx_range_x, [&](IdxX idx_x) {
            ddc::for_each(y_deriv_block, [&](Idx<dY> idx_dy) {
                dy_vals(idx_dy - y_deriv_block.front(),
                        idx_x - idx_range_x.front(),
                        deriv_idx_range_y.get_index(idx_y))
                        = example_grid(Idx_dXdYXY(idx_x, idx_y, no_deriv_idx_x, idx_dy));
            });
        });
    }

    // Get an mdspan describing the cross-derivatives (d_xd_y, d_xd_y^2, d_x^2d_y, d_x^2d_y^2) of the function
    detail::ViewNDMaker<4, int, false>::type dx_dy_vals = dxdyField.get_mdspan(x_y_deriv_block);
    // Check the shape of the mdspan
    EXPECT_EQ(dx_dy_vals.rank(), 4);
    EXPECT_EQ(dx_dy_vals.extent(0), x_deriv_block.size());
    EXPECT_EQ(dx_dy_vals.extent(1), y_deriv_block.size());
    EXPECT_EQ(dx_dy_vals.extent(2), deriv_idx_range_x.size());
    EXPECT_EQ(dx_dy_vals.extent(3), deriv_idx_range_y.size());

    // Fill the mdspan with values
    ddc::for_each(x_deriv_block, [&](Idx<dX> idx_dx) {
        ddc::for_each(y_deriv_block, [&](Idx<dY> idx_dy) {
            for (IdxX idx_x : deriv_idx_range_x) {
                for (IdxY idx_y : deriv_idx_range_y) {
                    dx_dy_vals(
                            idx_dx - x_deriv_block.front(),
                            idx_dy - y_deriv_block.front(),
                            deriv_idx_range_x.get_index(idx_x),
                            deriv_idx_range_y.get_index(idx_y))
                            = example_grid(Idx_dXdYXY(idx_x, idx_y, idx_dx, idx_dy));
                }
            }
        });
    });

    // Check that the value of the function at the index val_element is correct
    IdxXY val_element(lbound_x + 2, lbound_y + 4);
    EXPECT_EQ(dxdyField(val_element), example_grid(Idx_dXdYXY(val_element, no_deriv_idx)));

    // Check that the value of the 0-th derivatives of the function at the index val_element is correct
    Idx_dXdYXY exact_val_element(Idx<dX>(0), Idx<dY>(0), lbound_x + 2, lbound_y + 4);
    EXPECT_EQ(dxdyField(exact_val_element), example_grid(Idx_dXdYXY(exact_val_element)));

    // Check that the value of the 1st x-derivative of the function at the index (lbound_x, lbound_y) is correct
    Idx<dX, GridX, GridY> dx_element(Idx<dX>(1), lbound_x, lbound_y);
    EXPECT_EQ(dxdyField(dx_element), example_grid(Idx_dXdYXY(dx_element, no_deriv_idx_y)));

    // Check that the value of the 1st y-derivative of the function at the index (lbound_x, lbound_y) is correct
    Idx<dY, GridX, GridY> dy_element(Idx<dY>(1), lbound_x, lbound_y);
    EXPECT_EQ(dxdyField(dy_element), example_grid(Idx_dXdYXY(dy_element, no_deriv_idx_x)));

    // Check that the value of the cross derivative d_xd_y of the function at the index (lbound_x, lbound_y) is correct
    Idx_dXdYXY dx_dy_element(Idx<dX>(1), Idx<dY>(1), lbound_x, lbound_y);
    EXPECT_EQ(dxdyField(dx_dy_element), example_grid(dx_dy_element));
}

// Test the element-wise operator on GPU
void test_DerivField_GPUElementAccess()
{
    // Index ranges where derivatives are defined
    IdxRangeSlice<GridX> deriv_idx_range_x(idx_range_x.front(), IdxStepX(2), idx_range_x.extents());
    IdxRangeSlice<GridY> deriv_idx_range_y(idx_range_y.front(), IdxStepY(2), idx_range_y.extents());

    // Define a field memory allocation on x-y with 3 derivatives in x and y on GPU
    device_t<DerivFieldMem<int, IdxRange_dXdYXY, 3>>
            dxdyField_alloc(idx_range_x_y, deriv_idx_range_x, deriv_idx_range_y);
    // Define the field on GPU
    device_t<DerivField<int, IdxRange_dXdYXY>> dxdyField(dxdyField_alloc);

    // A subset  of the x-derivatives to be retrieved with get_mdspan
    IdxRange<dX> x_deriv_block(Idx<dX>(1), IdxStep<dX>(2));
    IdxRange<dY> y_deriv_block(Idx<dY>(1), IdxStep<dY>(2));
    IdxRange<dX, dY> x_y_deriv_block(x_deriv_block, y_deriv_block);

    // Get an mdspan describing the values of the function
    detail::ViewNDMaker<2, int, false>::type vals = dxdyField.get_mdspan();
    // Check the shape of the mdspan
    EXPECT_EQ(vals.rank(), 2);
    EXPECT_EQ(vals.extent(0), idx_range_x.size());
    EXPECT_EQ(vals.extent(1), idx_range_y.size());

    GridBuilder example_grid(idx_range_x_y_dx_dy);

    Idx<dX, dY> no_deriv_idx(0, 0);
    // Fill the mdspan with values on GPU
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            idx_range_x_y,
            KOKKOS_LAMBDA(IdxXY idx_x_y) {
                Idx<GridX> idx_x(idx_x_y);
                Idx<GridY> idx_y(idx_x_y);
                vals(idx_x - idx_range_x.front(), idx_y - idx_range_y.front())
                        = example_grid(Idx_dXdYXY(idx_x_y, no_deriv_idx));
            });

    // Get an mdspan describing the 1st and 2nd x-derivatives of the function
    detail::ViewNDMaker<3, int, false>::type dx_vals = dxdyField.get_mdspan(x_deriv_block);
    // Check the shape of the mdspan
    EXPECT_EQ(dx_vals.rank(), 3);
    EXPECT_EQ(dx_vals.extent(0), x_deriv_block.size());
    EXPECT_EQ(dx_vals.extent(1), deriv_idx_range_x.size());
    EXPECT_EQ(dx_vals.extent(2), idx_range_y.size());

    // Collapse the index ranges into 1 iterable
    IdxRange<dX, GridY> idx_range_dx_y(x_deriv_block, idx_range_y);

    Idx<dY> no_deriv_idx_y(0);
    // Fill the mdspan with values on GPU
    for (IdxX idx_x : deriv_idx_range_x) {
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                idx_range_dx_y,
                KOKKOS_LAMBDA(Idx<dX, GridY> idx_dx_y) {
                    Idx<dX> idx_dx(idx_dx_y);
                    Idx<GridY> idx_y(idx_dx_y);
                    dx_vals(idx_dx - x_deriv_block.front(),
                            deriv_idx_range_x.get_index(idx_x),
                            idx_y - idx_range_y.front())
                            = example_grid(Idx_dXdYXY(idx_x, idx_y, idx_dx, no_deriv_idx_y));
                });
    }

    // Get an mdspan describing the 1st and 2nd y-derivatives of the function
    detail::ViewNDMaker<3, int, false>::type dy_vals = dxdyField.get_mdspan(y_deriv_block);
    // Check the shape of the mdspan
    EXPECT_EQ(dy_vals.rank(), 3);
    EXPECT_EQ(dy_vals.extent(0), y_deriv_block.size());
    EXPECT_EQ(dy_vals.extent(1), idx_range_x.size());
    EXPECT_EQ(dy_vals.extent(2), deriv_idx_range_y.size());

    // Collapse the index ranges into 1 iterable
    IdxRange<dY, GridX> idx_range_x_dy(y_deriv_block, idx_range_x);

    Idx<dX> no_deriv_idx_x(0);
    // Fill the mdspan with values on GPU
    for (IdxY idx_y : deriv_idx_range_y) {
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                idx_range_x_dy,
                KOKKOS_LAMBDA(Idx<dY, GridX> idx_x_dy) {
                    Idx<dY> idx_dy(idx_x_dy);
                    Idx<GridX> idx_x(idx_x_dy);
                    dy_vals(idx_dy - y_deriv_block.front(),
                            idx_x - idx_range_x.front(),
                            deriv_idx_range_y.get_index(idx_y))
                            = example_grid(Idx_dXdYXY(idx_x, idx_y, no_deriv_idx_x, idx_dy));
                });
    }

    // Get an mdspan describing the cross-derivatives (d_xd_y, d_xd_y^2, d_x^2d_y, d_x^2d_y) of the function
    detail::ViewNDMaker<4, int, false>::type dx_dy_vals = dxdyField.get_mdspan(x_y_deriv_block);
    // Check the shape of the mdspan
    EXPECT_EQ(dx_dy_vals.rank(), 4);
    EXPECT_EQ(dx_dy_vals.extent(0), x_deriv_block.size());
    EXPECT_EQ(dx_dy_vals.extent(1), y_deriv_block.size());
    EXPECT_EQ(dx_dy_vals.extent(2), deriv_idx_range_x.size());
    EXPECT_EQ(dx_dy_vals.extent(3), deriv_idx_range_y.size());

    // Collapse the index ranges into 1 iterable
    IdxRange<dX, dY> idx_range_dx_dy(x_deriv_block, y_deriv_block);

    // Fill the mdspan with values on GPU
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            idx_range_dx_dy,
            KOKKOS_LAMBDA(Idx<dX, dY> idx_dx_dy) {
                Idx<dX> idx_dx(idx_dx_dy);
                Idx<dY> idx_dy(idx_dx_dy);
                for (IdxX idx_x : deriv_idx_range_x) {
                    for (IdxY idx_y : deriv_idx_range_y) {
                        dx_dy_vals(
                                idx_dx - x_deriv_block.front(),
                                idx_dy - y_deriv_block.front(),
                                deriv_idx_range_x.get_index(idx_x),
                                deriv_idx_range_y.get_index(idx_y))
                                = example_grid(Idx_dXdYXY(idx_x, idx_y, idx_dx, idx_dy));
                    }
                }
            });

    // Copy to non-strided layout to be able to call create_mirror_view_and_copy
    FieldMem<int, IdxRangeXY> ddc_vals_alloc(idx_range_x_y);
    Field<int, IdxRangeXY> ddc_vals = get_field(ddc_vals_alloc);
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            idx_range_x_y,
            KOKKOS_LAMBDA(IdxXY idx_x_y) { ddc_vals(idx_x_y) = dxdyField(idx_x_y); });

    // Check values
    auto vals_host = ddc::create_mirror_view_and_copy(get_const_field(ddc_vals));
    ddc::for_each(idx_range_x_y, [&](IdxXY idx_x_y) {
        EXPECT_EQ(vals_host(idx_x_y), example_grid(Idx_dXdYXY(idx_x_y, no_deriv_idx)));
    });

    IdxRange<GridX> idx_range_x_slice(idx_range_x.front(), IdxStep<GridX>(1));
    IdxRange<dX, GridX, GridY>
            idx_range_dx_x_y_slice(x_deriv_block, idx_range_x_slice, idx_range_y);
    // Copy to non-strided layout to be able to call create_mirror_view_and_copy
    FieldMem<int, IdxRange<dX, GridX, GridY>> ddc_dx_vals_alloc(idx_range_dx_x_y_slice);
    Field<int, IdxRange<dX, GridX, GridY>> ddc_dx_vals = get_field(ddc_dx_vals_alloc);
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            idx_range_dx_x_y_slice,
            KOKKOS_LAMBDA(Idx<dX, GridX, GridY> idx) { ddc_dx_vals(idx) = dxdyField(idx); });

    // Check x-derivative values
    auto dx_vals_host = ddc::create_mirror_view_and_copy(get_const_field(ddc_dx_vals));
    ddc::for_each(idx_range_dx_x_y_slice, [&](Idx<dX, GridX, GridY> idx_dx_x_y) {
        Idx<dX> idx_dx(idx_dx_x_y);
        Idx<GridX> idx_x(idx_dx_x_y);
        Idx<GridY> idx_y(idx_dx_x_y);
        EXPECT_EQ(dx_vals_host(idx_dx_x_y), example_grid(Idx_dXdYXY(idx_dx_x_y, no_deriv_idx_y)));
    });

    IdxRange<GridY> idx_range_y_slice(idx_range_y.front(), IdxStep<GridY>(1));
    IdxRange<dY, GridX, GridY>
            idx_range_dy_x_y_slice(y_deriv_block, idx_range_x, idx_range_y_slice);
    // Copy to non-strided layout to be able to call create_mirror_view_and_copy
    FieldMem<int, IdxRange<dY, GridX, GridY>> ddc_dy_vals_alloc(idx_range_dy_x_y_slice);
    Field<int, IdxRange<dY, GridX, GridY>> ddc_dy_vals = get_field(ddc_dy_vals_alloc);
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            idx_range_dy_x_y_slice,
            KOKKOS_LAMBDA(Idx<dY, GridX, GridY> idx) { ddc_dy_vals(idx) = dxdyField(idx); });

    // Check y-derivative values
    auto dy_vals_host = ddc::create_mirror_view_and_copy(get_const_field(ddc_dy_vals));
    ddc::for_each(idx_range_dy_x_y_slice, [&](Idx<dY, GridX, GridY> idx_dy_x_y) {
        Idx<dY> idx_dy(idx_dy_x_y);
        Idx<GridX> idx_x(idx_dy_x_y);
        Idx<GridY> idx_y(idx_dy_x_y);
        EXPECT_EQ(dy_vals_host(idx_dy_x_y), example_grid(Idx_dXdYXY(idx_dy_x_y, no_deriv_idx_x)));
    });

    IdxRange_dXdYXY idx_range_dx_dy_x_y_slice(
            x_deriv_block,
            y_deriv_block,
            idx_range_x_slice,
            idx_range_y_slice);
    // Copy to non-strided layout to be able to call create_mirror_view_and_copy
    FieldMem<int, IdxRange_dXdYXY> ddc_dx_dy_vals_alloc(idx_range_dx_dy_x_y_slice);
    Field<int, IdxRange_dXdYXY> ddc_dx_dy_vals = get_field(ddc_dx_dy_vals_alloc);
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            idx_range_dx_dy_x_y_slice,
            KOKKOS_LAMBDA(Idx_dXdYXY idx) { ddc_dx_dy_vals(idx) = dxdyField(idx); });

    // Check cross-derivative values
    auto dx_dy_vals_host = ddc::create_mirror_view_and_copy(get_const_field(ddc_dx_dy_vals));
    ddc::for_each(idx_range_dx_dy_x_y_slice, [&](Idx_dXdYXY idx_dx_dy_x_y) {
        EXPECT_EQ(dx_dy_vals_host(idx_dx_dy_x_y), example_grid(Idx_dXdYXY(idx_dx_dy_x_y)));
    });
}

TEST(DerivFieldTest, GPUElementAccess)
{
    test_DerivField_GPUElementAccess();
}

// Test if the deepcopy correctly fills a field mem
TEST(DerivFieldMemTest, FieldDeepCopy)
{
    // Index ranges where derivatives are defined
    IdxRangeSlice<GridX> deriv_idx_range_x(idx_range_x.front(), IdxStepX(2), idx_range_x.extents());

    // Define fields on x-y with 1 derivative in x and y
    DerivFieldMem<double, IdxRange<dX, GridX, GridY>, 1>
            dxdyField(idx_range_x_y, deriv_idx_range_x);
    DerivFieldMem<double, IdxRange<dX, GridX, GridY>, 1>
            dxdyField_copy(idx_range_x_y, deriv_idx_range_x);

    // Extract the values and derivatives
    Idx<dX> first_deriv(1);
    host_t<DField<IdxRange<GridX, GridY>, std::experimental::layout_stride>> values
            = dxdyField.get_values_field();
    host_t<DField<IdxRange<GridY>, std::experimental::layout_stride>> left_x_derivs
            = dxdyField[IdxXdX(first_deriv, deriv_idx_range_x.front())];
    host_t<DField<IdxRange<GridY>, std::experimental::layout_stride>> right_x_derivs
            = dxdyField[IdxXdX(first_deriv, deriv_idx_range_x.back())];

    // Set the values
    ddc::for_each(get_idx_range(values), [&](IdxXY ixy) {
        IdxX ix(ixy);
        IdxY iy(ixy);
        double x(ix - IdxX(0));
        double y(0.1 * (iy - IdxY(0)));
        values(ixy) = x * x + y * y;
    });
    // Set the derivatives
    double x_left(idx_range_x.front() - IdxX(0));
    ddc::for_each(get_idx_range(left_x_derivs), [&](IdxY iy) { left_x_derivs(iy) = 2 * x_left; });
    double x_right(idx_range_x.back() - IdxX(0));
    ddc::for_each(get_idx_range(right_x_derivs), [&](IdxY iy) {
        right_x_derivs(iy) = 2 * x_right;
    });

    // Copy the deriv field
    ddcHelper::deepcopy(dxdyField_copy, dxdyField);

    // Extract the values and derivatives from the copy
    host_t<DField<IdxRange<GridX, GridY>, std::experimental::layout_stride>> values_copy
            = dxdyField_copy.get_values_field();
    host_t<DField<IdxRange<GridY>, std::experimental::layout_stride>> left_x_derivs_copy
            = dxdyField[IdxXdX(first_deriv, deriv_idx_range_x.front())];
    host_t<DField<IdxRange<GridY>, std::experimental::layout_stride>> right_x_derivs_copy
            = dxdyField[IdxXdX(first_deriv, deriv_idx_range_x.back())];

    // Check the values
    ddc::for_each(get_idx_range(values), [&](IdxXY ixy) {
        IdxX ix(ixy);
        IdxY iy(ixy);
        double x(ix - IdxX(0));
        double y(0.1 * (iy - IdxY(0)));
        EXPECT_EQ(values_copy(ixy), x * x + y * y);
    });
    // Check the derivatives
    ddc::for_each(get_idx_range(left_x_derivs), [&](IdxY iy) {
        EXPECT_EQ(left_x_derivs(iy), 2 * x_left);
    });
    ddc::for_each(get_idx_range(right_x_derivs), [&](IdxY iy) {
        EXPECT_EQ(right_x_derivs(iy), 2 * x_right);
    });
}

// Test if the deepcopy correctly fills a field
TEST(DerivFieldTest, FieldDeepCopy)
{
    // Index ranges where derivatives are defined
    IdxRangeSlice<GridX> deriv_idx_range_x(idx_range_x.front(), IdxStepX(2), idx_range_x.extents());

    // Define field memory allocations on x-y with 1 derivative in x and y
    DerivFieldMem<double, IdxRange<dX, GridX, GridY>, 1>
            dxdyField_alloc(idx_range_x_y, deriv_idx_range_x);
    DerivFieldMem<double, IdxRange<dX, GridX, GridY>, 1>
            dxdyField_copy_alloc(idx_range_x_y, deriv_idx_range_x);

    // Get the fields
    DerivField dxdyField = get_field(dxdyField_alloc);
    DerivField dxdyField_copy = get_field(dxdyField_copy_alloc);

    // Extract the values and derivatives
    Idx<dX> first_deriv(1);
    host_t<DField<IdxRange<GridX, GridY>, std::experimental::layout_stride>> values
            = dxdyField.get_values_field();
    host_t<DField<IdxRange<GridY>, std::experimental::layout_stride>> left_x_derivs
            = dxdyField[IdxXdX(first_deriv, deriv_idx_range_x.front())];
    host_t<DField<IdxRange<GridY>, std::experimental::layout_stride>> right_x_derivs
            = dxdyField[IdxXdX(first_deriv, deriv_idx_range_x.back())];

    // Set the values
    ddc::for_each(get_idx_range(values), [&](IdxXY ixy) {
        IdxX ix(ixy);
        IdxY iy(ixy);
        double x(ix - IdxX(0));
        double y(0.1 * (iy - IdxY(0)));
        values(ixy) = x * x + y * y;
    });
    // Set the derivatives
    double x_left(idx_range_x.front() - IdxX(0));
    ddc::for_each(get_idx_range(left_x_derivs), [&](IdxY iy) { left_x_derivs(iy) = 2 * x_left; });
    double x_right(idx_range_x.back() - IdxX(0));
    ddc::for_each(get_idx_range(right_x_derivs), [&](IdxY iy) {
        right_x_derivs(iy) = 2 * x_right;
    });

    // Copy the deriv field
    ddcHelper::deepcopy(dxdyField_copy, dxdyField);

    // Extract the values and derivatives from the copy
    host_t<DField<IdxRange<GridX, GridY>, std::experimental::layout_stride>> values_copy
            = dxdyField_copy.get_values_field();
    host_t<DField<IdxRange<GridY>, std::experimental::layout_stride>> left_x_derivs_copy
            = dxdyField[IdxXdX(first_deriv, deriv_idx_range_x.front())];
    host_t<DField<IdxRange<GridY>, std::experimental::layout_stride>> right_x_derivs_copy
            = dxdyField[IdxXdX(first_deriv, deriv_idx_range_x.back())];

    // Check the values
    ddc::for_each(get_idx_range(values), [&](IdxXY ixy) {
        IdxX ix(ixy);
        IdxY iy(ixy);
        double x(ix - IdxX(0));
        double y(0.1 * (iy - IdxY(0)));
        EXPECT_EQ(values(ixy), x * x + y * y);
    });
    // Check the derivatives
    ddc::for_each(get_idx_range(left_x_derivs), [&](IdxY iy) {
        EXPECT_EQ(left_x_derivs(iy), 2 * x_left);
    });
    ddc::for_each(get_idx_range(right_x_derivs), [&](IdxY iy) {
        EXPECT_EQ(right_x_derivs(iy), 2 * x_right);
    });
}
