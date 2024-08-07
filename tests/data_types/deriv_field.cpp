// SPDX-License-Identifier: MIT
#include <ddc/ddc.hpp>
#include <ddc/kernels/splines/deriv.hpp>

#include <gtest/gtest.h>

#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "derivative_field.hpp"
#include "derivative_field_span.hpp"
#include "directional_tag.hpp"

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

static IdxX constexpr lbound_x(50);
static IdxStepX constexpr nelems_x(3);
static IdxRangeX constexpr dom_x(lbound_x, nelems_x);

static IdxY constexpr lbound_y(4);
static IdxStepY constexpr nelems_y(12);
static IdxRangeY constexpr dom_y(lbound_y, nelems_y);

static IdxXY constexpr lbound_x_y {lbound_x, lbound_y};
static IdxStepXY constexpr nelems_x_y(nelems_x, nelems_y);
static IdxRangeXY constexpr dom_x_y(lbound_x_y, nelems_x_y);
} // namespace

// Test the constructor for a 2D field with derivatives in 1 direction
TEST(DerivFieldTest, Constructor1Deriv)
{
    // Type for a x,y field with 1 derivative in x
    using DFieldXY_dX = DerivField<double, IdxRange<dX, GridX, GridY>, 1>;

    // Domain where derivatives are defined
    DiscreteSubDomain<GridX> deriv_idx_range_x(dom_x.front(), IdxStepX(2), dom_x.extents());

    // Define the field
    DFieldXY_dX dxField(dom_x_y, deriv_idx_range_x);

    // Ensure that the internal chunk has the expected type
    constexpr bool same = std::is_same_v<
            typename DFieldXY_dX::chunk_type,
            host_t<DFieldMem<IdxRange<dX, GridX, GridY>>>>;
    EXPECT_TRUE(same);
}

// Test the constructor for a 2D field with derivatives in 2 directions
TEST(DerivFieldTest, Constructor2Deriv)
{
    // Type for a x,y field with 1 derivative in x and 1 derivative in y
    using DFieldXY_dXdY = DerivField<double, IdxRange<dX, dY, GridX, GridY>, 1>;

    // Domains where derivatives are defined
    DiscreteSubDomain<GridX> deriv_idx_range_x(dom_x.front(), IdxStepX(2), dom_x.extents());
    DiscreteSubDomain<GridY> deriv_idx_range_y(dom_y.front(), IdxStepY(2), dom_y.extents());

    // Define the field
    DFieldXY_dXdY dxdyField(dom_x_y, deriv_idx_range_x, deriv_idx_range_y);

    // Ensure that the internal chunk has the expected type
    bool same = std::is_same_v<
            typename DFieldXY_dXdY::chunk_type,
            host_t<DFieldMem<IdxRange<dX, dY, GridX, GridY>>>>;
    EXPECT_TRUE(same);
}

// Test the span constructor for a 2D field with derivatives in 1 direction
TEST(DerivFieldSpanTest, Constructor1Deriv)
{
    // Type for a x,y field span with 1 derivative in x
    using DSpanXY_dX = DerivFieldSpan<double, IdxRange<dX, GridX, GridY>>;

    // Domain where derivatives are defined
    DiscreteSubDomain<GridX> deriv_idx_range_x(dom_x.front(), IdxStepX(2), dom_x.extents());

    // Define a field on x-y with 1 derivative in x and y
    DerivField<double, IdxRange<dX, GridX, GridY>, 1> dxField_alloc(dom_x_y, deriv_idx_range_x);
    // Define the span
    DSpanXY_dX dxField(dxField_alloc);

    // Ensure that the internal chunk span has the expected type
    bool same = std::is_same_v<DSpanXY_dX::chunk_type, host_t<DField<IdxRange<dX, GridX, GridY>>>>;
    EXPECT_TRUE(same);
}

// Test the span constructor for a 2D field with derivatives in 2 directions
TEST(DerivFieldSpanTest, Constructor2Deriv)
{
    // Type for a x,y field span with 1 derivative in x and 1 derivative in y
    using DSpanXY_dXdY = DerivFieldSpan<double, IdxRange<dX, dY, GridX, GridY>>;

    // Domains where derivatives are defined
    DiscreteSubDomain<GridX> deriv_idx_range_x(dom_x.front(), IdxStepX(2), dom_x.extents());
    DiscreteSubDomain<GridY> deriv_idx_range_y(dom_y.front(), IdxStepY(2), dom_y.extents());

    // Define a field on x-y with 1 derivative in x and y
    DerivField<double, IdxRange<dX, dY, GridX, GridY>, 1>
            dxdyField_alloc(dom_x_y, deriv_idx_range_x, deriv_idx_range_y);
    // Define the span
    DSpanXY_dXdY dxdyField(dxdyField_alloc);

    // Ensure that the internal chunk span has the expected type
    bool same = std::
            is_same_v<DSpanXY_dXdY::chunk_type, host_t<DField<IdxRange<dX, dY, GridX, GridY>>>>;
    EXPECT_TRUE(same);
}

// Test if the values of the function can be accessed via the get_values_span function
TEST(DerivFieldSpanTest, SpanValueAccess)
{
    // Domains where derivatives are defined
    DiscreteSubDomain<GridX> deriv_idx_range_x(dom_x.front(), IdxStepX(2), dom_x.extents());
    DiscreteSubDomain<GridY> deriv_idx_range_y(dom_y.front(), IdxStepY(2), dom_y.extents());

    // Define a field on x-y with 1 derivative in x and y
    DerivField<double, IdxRange<dX, dY, GridX, GridY>, 1>
            dxdyField_alloc(dom_x_y, deriv_idx_range_x, deriv_idx_range_y);
    // Define the span
    DerivFieldSpan<double, IdxRange<dX, dY, GridX, GridY>> dxdyField(dxdyField_alloc);

    // Check that the index range of the values matches the expected index range
    EXPECT_EQ(dom_x_y, get_idx_range(dxdyField.get_values_span()));
}

// Test if the values of the function can be accessed via the constant get_values_span function
TEST(DerivFieldSpanTest, ViewValueAccess)
{
    // Domains where derivatives are defined
    DiscreteSubDomain<GridX> deriv_idx_range_x(dom_x.front(), IdxStepX(2), dom_x.extents());
    DiscreteSubDomain<GridY> deriv_idx_range_y(dom_y.front(), IdxStepY(2), dom_y.extents());

    // Define a field on x-y with 1 derivative in x and y
    DerivField<double, IdxRange<dX, dY, GridX, GridY>, 1>
            dxdyField_alloc(dom_x_y, deriv_idx_range_x, deriv_idx_range_y);
    // Define the span
    DerivFieldView<const double, IdxRange<dX, dY, GridX, GridY>> dxdyField(dxdyField_alloc);

    // Check that the index range of the values matches the expected index range
    EXPECT_EQ(dom_x_y, get_idx_range(dxdyField.get_values_span()));
}


// Test if the derivatives of the function can be accessed via the slice function
TEST(DerivFieldTest, derivValueAccess)
{
    // Domains where derivatives are defined
    DiscreteSubDomain<GridX> deriv_idx_range_x(dom_x.front(), IdxStepX(2), dom_x.extents());
    DiscreteSubDomain<GridY> deriv_idx_range_y(dom_y.front(), IdxStepY(2), dom_y.extents());

    // Define a field on x-y with 2 derivatives in x and y
    DerivField<double, IdxRange<dX, dY, GridX, GridY>, 2>
            dxdyField(dom_x_y, deriv_idx_range_x, deriv_idx_range_y);

    // The derivatives to be retrieved from the field
    IdxRange<dX> x_deriv_block(Idx<dX>(1), IdxStep<dX>(2));
    IdxRange<dY> y_deriv_block(Idx<dY>(1), IdxStep<dY>(2));

    // The subindex ranges to be retrieved from the field
    IdxRange<GridX> deriv_x_subdom(dom_x.front(), IdxStep<GridX>(1));
    IdxRange<GridY> deriv_y_subdom(dom_y.front(), IdxStep<GridY>(1));

    // The expected index range for the x-derivatives identified in x_deriv_block at the positions identified by deriv_x_subdom
    IdxRange<dX, GridX, GridY>
            dx_idx_range(x_deriv_block, deriv_x_subdom, ddc::select<GridY>(dom_x_y));
    // The expected index range for the y-derivatives identified in y_deriv_block at the positions identified by deriv_y_subdom
    IdxRange<dY, GridX, GridY>
            dy_idx_range(y_deriv_block, ddc::select<GridX>(dom_x_y), deriv_y_subdom);
    // The expected index range for the x-y-derivatives identified in x_deriv_block and y_deriv_block at the positions identified
    // by deriv_x_subdom and deriv_y_subdom
    IdxRange<dX, dY, GridX, GridY>
            dx_dy_idx_range(x_deriv_block, y_deriv_block, deriv_x_subdom, deriv_y_subdom);

    // Build the multi-D index ranges that will slice the field
    IdxRange<dX, GridX> slice_idx_dx(x_deriv_block, deriv_x_subdom);
    IdxRange<dY, GridY> slice_idx_dy(y_deriv_block, deriv_y_subdom);
    IdxRange<dX, dY, GridX, GridY>
            slice_idx_dx_dy(x_deriv_block, y_deriv_block, deriv_x_subdom, deriv_y_subdom);

    // Collect the index range of the sliced fields
    IdxRange<dX, GridX, GridY> slice_dx_idx_range = get_idx_range(dxdyField[slice_idx_dx]);
    IdxRange<dY, GridX, GridY> slice_dy_idx_range = get_idx_range(dxdyField[slice_idx_dy]);
    IdxRange<dX, dY, GridX, GridY> slice_dx_dy_idx_range
            = get_idx_range(dxdyField[slice_idx_dx_dy]);

    // Check that the index ranges are as expected
    EXPECT_EQ(dx_idx_range, slice_dx_idx_range);
    EXPECT_EQ(dy_idx_range, slice_dy_idx_range);
    EXPECT_EQ(dx_dy_idx_range, slice_dx_dy_idx_range);
}

// Test if the derivatives of the function can be accessed from a span via the slice function
TEST(DerivFieldSpanTest, derivSpanValueAccess)
{
    // Domains where derivatives are defined
    DiscreteSubDomain<GridX> deriv_idx_range_x(dom_x.front(), IdxStepX(2), dom_x.extents());
    DiscreteSubDomain<GridY> deriv_idx_range_y(dom_y.front(), IdxStepY(2), dom_y.extents());

    // Define a field on x-y with 2 derivatives in x and y
    DerivField<double, IdxRange<dX, dY, GridX, GridY>, 2>
            dxdyField_alloc(dom_x_y, deriv_idx_range_x, deriv_idx_range_y);
    // Define the span
    DerivFieldSpan<double, IdxRange<dX, dY, GridX, GridY>> dxdyField(dxdyField_alloc);

    // The derivatives to be retrieved from the field
    IdxRange<dX> x_deriv_block(Idx<dX>(1), IdxStep<dX>(2));
    IdxRange<dY> y_deriv_block(Idx<dY>(1), IdxStep<dY>(2));

    // The subindex ranges to be retrieved from the field
    IdxRange<GridX> deriv_x_subdom(dom_x.front(), IdxStep<GridX>(1));
    IdxRange<GridY> deriv_y_subdom(dom_y.front(), IdxStep<GridY>(1));

    // The expected index range for the x-derivatives identified in x_deriv_block at the positions identified by deriv_x_subdom
    IdxRange<dX, GridX, GridY>
            dx_idx_range(x_deriv_block, deriv_x_subdom, ddc::select<GridY>(dom_x_y));
    // The expected index range for the y-derivatives identified in y_deriv_block at the positions identified by deriv_y_subdom
    IdxRange<dY, GridX, GridY>
            dy_idx_range(y_deriv_block, ddc::select<GridX>(dom_x_y), deriv_y_subdom);
    // The expected index range for the x-y-derivatives identified in x_deriv_block and y_deriv_block at the positions identified
    // by deriv_x_subdom and deriv_y_subdom
    IdxRange<dX, dY, GridX, GridY>
            dx_dy_idx_range(x_deriv_block, y_deriv_block, deriv_x_subdom, deriv_y_subdom);

    // Build the multi-D index ranges that will slice the field
    IdxRange<dX, GridX> slice_idx_dx(x_deriv_block, deriv_x_subdom);
    IdxRange<dY, GridY> slice_idx_dy(y_deriv_block, deriv_y_subdom);
    IdxRange<dX, dY, GridX, GridY>
            slice_idx_dx_dy(x_deriv_block, y_deriv_block, deriv_x_subdom, deriv_y_subdom);

    // Collect the index range of the sliced fields
    IdxRange<dX, GridX, GridY> slice_dx_idx_range = get_idx_range(dxdyField[slice_idx_dx]);
    IdxRange<dY, GridX, GridY> slice_dy_idx_range = get_idx_range(dxdyField[slice_idx_dy]);
    IdxRange<dX, dY, GridX, GridY> slice_dx_dy_idx_range
            = get_idx_range(dxdyField[slice_idx_dx_dy]);

    // Check that the index ranges are as expected
    EXPECT_EQ(dx_idx_range, slice_dx_idx_range);
    EXPECT_EQ(dy_idx_range, slice_dy_idx_range);
    EXPECT_EQ(dx_dy_idx_range, slice_dx_dy_idx_range);
}

// Test the element-wise operator
TEST(DerivFieldSpanTest, ElementAccess)
{
    // Domains where derivatives are defined
    DiscreteSubDomain<GridX> deriv_idx_range_x(dom_x.front(), IdxStepX(2), dom_x.extents());
    DiscreteSubDomain<GridY> deriv_idx_range_y(dom_y.front(), IdxStepY(2), dom_y.extents());

    // Define a field on x-y with 3 derivatives in x and y
    DerivField<int, IdxRange<dX, dY, GridX, GridY>, 3>
            dxdyField_alloc(dom_x_y, deriv_idx_range_x, deriv_idx_range_y);
    // Define the span
    DerivFieldSpan<int, IdxRange<dX, dY, GridX, GridY>> dxdyField(dxdyField_alloc);

    // A subset  of the x-derivatives to be retrieved with get_mdspan
    IdxRange<dX> x_deriv_block(Idx<dX>(1), IdxStep<dX>(2));
    // A subset  of the y-derivatives to be retrieved with get_mdspan
    IdxRange<dY> y_deriv_block(Idx<dY>(1), IdxStep<dY>(2));
    // A subset  of the x and y derivatives to be retrieved with get_mdspan
    IdxRange<dX, dY> x_y_deriv_block(x_deriv_block, y_deriv_block);

    // Get an mdspan describing the values of the function
    auto vals = dxdyField.get_mdspan();
    // Check the shape of the span
    EXPECT_EQ(vals.rank(), 2);
    EXPECT_EQ(vals.extent(0), dom_x.size());
    EXPECT_EQ(vals.extent(1), dom_y.size());

    // Fill the span with values
    ddc::for_each(dom_x, [&](auto idx_x) {
        ddc::for_each(dom_y, [&](auto idx_y) {
            vals(idx_x - dom_x.front(), idx_y - dom_y.front())
                    = idx_x.uid() * nelems_y.value() + idx_y.uid();
        });
    });

    // Get an mdspan describing the 1st and 2nd x-derivatives of the function
    auto dx_vals = dxdyField.get_mdspan(x_deriv_block);
    // Check the shape of the span
    EXPECT_EQ(dx_vals.rank(), 3);
    EXPECT_EQ(dx_vals.extent(0), x_deriv_block.size());
    EXPECT_EQ(dx_vals.extent(1), deriv_idx_range_x.size());
    EXPECT_EQ(dx_vals.extent(2), dom_y.size());

    // Fill the span with values
    for (auto idx_x : deriv_idx_range_x) {
        ddc::for_each(dom_y, [&](auto idx_y) {
            ddc::for_each(x_deriv_block, [&](auto idx_dx) {
                dx_vals(idx_dx - x_deriv_block.front(),
                        deriv_idx_range_x.get_index(idx_x),
                        idx_y - dom_y.front())
                        = 100 + idx_dx.uid() * 50 + idx_x.uid() * nelems_y.value() + idx_y.uid();
            });
        });
    }

    // Get an mdspan describing the 1st and 2nd y-derivatives of the function
    auto dy_vals = dxdyField.get_mdspan(y_deriv_block);
    // Check the shape of the span
    EXPECT_EQ(dy_vals.rank(), 3);
    EXPECT_EQ(dy_vals.extent(0), y_deriv_block.size());
    EXPECT_EQ(dy_vals.extent(1), dom_x.size());
    EXPECT_EQ(dy_vals.extent(2), deriv_idx_range_y.size());

    // Fill the span with values
    for (auto idx_y : deriv_idx_range_y) {
        ddc::for_each(dom_x, [&](auto idx_x) {
            ddc::for_each(y_deriv_block, [&](auto idx_dy) {
                dy_vals(idx_dy - y_deriv_block.front(),
                        idx_x - dom_x.front(),
                        deriv_idx_range_y.get_index(idx_y))
                        = 200 + idx_dy.uid() * 50 + idx_x.uid() * nelems_y.value() + idx_y.uid();
            });
        });
    }

    // Get an mdspan describing the cross-derivatives (d_xd_y, d_xd_y^2, d_x^2d_y, d_x^2d_y^2) of the function
    auto dx_dy_vals = dxdyField.get_mdspan(x_y_deriv_block);
    // Check the shape of the span
    EXPECT_EQ(dx_dy_vals.rank(), 4);
    EXPECT_EQ(dx_dy_vals.extent(0), x_deriv_block.size());
    EXPECT_EQ(dx_dy_vals.extent(1), y_deriv_block.size());
    EXPECT_EQ(dx_dy_vals.extent(2), deriv_idx_range_x.size());
    EXPECT_EQ(dx_dy_vals.extent(3), deriv_idx_range_y.size());

    // Fill the span with values
    ddc::for_each(x_deriv_block, [&](auto idx_dx) {
        ddc::for_each(y_deriv_block, [&](auto idx_dy) {
            for (auto idx_x : deriv_idx_range_x) {
                for (auto idx_y : deriv_idx_range_y) {
                    dx_dy_vals(
                            idx_dx - x_deriv_block.front(),
                            idx_dy - y_deriv_block.front(),
                            deriv_idx_range_x.get_index(idx_x),
                            deriv_idx_range_y.get_index(idx_y))
                            = 300 + idx_dx.uid() * 100 + idx_dy.uid() * 50
                              + idx_x.uid() * nelems_y.value() + idx_y.uid();
                }
            }
        });
    });

    // Check that the value of the function at the index val_element is correct
    IdxXY val_element(lbound_x + 2, lbound_y + 4);
    EXPECT_EQ(
            dxdyField(val_element),
            ddc::uid<GridX>(val_element) * nelems_y.value() + ddc::uid<GridY>(val_element));

    // Check that the value of the 0-th derivatives of the function at the index val_element is correct
    Idx<dX, dY, GridX, GridY> exact_val_element(Idx<dX>(0), Idx<dY>(0), lbound_x + 2, lbound_y + 4);
    EXPECT_EQ(
            dxdyField(exact_val_element),
            ddc::uid<GridX>(val_element) * nelems_y.value() + ddc::uid<GridY>(val_element));

    // Check that the value of the 1st x-derivative of the function at the index (lbound_x, lbound_y) is correct
    Idx<dX, GridX, GridY> dx_element(Idx<dX>(1), lbound_x, lbound_y);
    EXPECT_EQ(
            dxdyField(dx_element),
            100 + ddc::uid<dX>(dx_element) * 50 + ddc::uid<GridX>(dx_element) * nelems_y.value()
                    + ddc::uid<GridY>(dx_element));

    // Check that the value of the 1st y-derivative of the function at the index (lbound_x, lbound_y) is correct
    Idx<dY, GridX, GridY> dy_element(Idx<dY>(1), lbound_x, lbound_y);
    EXPECT_EQ(
            dxdyField(dy_element),
            200 + ddc::uid<dY>(dy_element) * 50 + ddc::uid<GridX>(dy_element) * nelems_y.value()
                    + ddc::uid<GridY>(dy_element));

    // Check that the value of the cross derivative d_xd_y of the function at the index (lbound_x, lbound_y) is correct
    Idx<dX, dY, GridX, GridY> dx_dy_element(Idx<dX>(1), Idx<dY>(1), lbound_x, lbound_y);
    EXPECT_EQ(
            dxdyField(dx_dy_element),
            300 + ddc::uid<dX>(dx_dy_element) * 100 + ddc::uid<dY>(dx_dy_element) * 50
                    + ddc::uid<GridX>(dx_dy_element) * nelems_y.value()
                    + ddc::uid<GridY>(dx_dy_element));
}

// Test the element-wise operator on GPU
void test_DerivField_GPUElementAccess()
{
    // Domains where derivatives are defined
    DiscreteSubDomain<GridX> deriv_idx_range_x(dom_x.front(), IdxStepX(2), dom_x.extents());
    DiscreteSubDomain<GridY> deriv_idx_range_y(dom_y.front(), IdxStepY(2), dom_y.extents());

    // Define a field on x-y with 3 derivatives in x and y on GPU
    device_t<DerivField<int, IdxRange<dX, dY, GridX, GridY>, 3>>
            dxdyField_alloc(dom_x_y, deriv_idx_range_x, deriv_idx_range_y);
    // Define the span on GPU
    device_t<DerivFieldSpan<int, IdxRange<dX, dY, GridX, GridY>>> dxdyField(dxdyField_alloc);

    // A subset  of the x-derivatives to be retrieved with get_mdspan
    IdxRange<dX> x_deriv_block(Idx<dX>(1), IdxStep<dX>(2));
    IdxRange<dY> y_deriv_block(Idx<dY>(1), IdxStep<dY>(2));
    IdxRange<dX, dY> x_y_deriv_block(x_deriv_block, y_deriv_block);

    // Get an mdspan describing the values of the function
    auto vals = dxdyField.get_mdspan();
    // Check the shape of the span
    EXPECT_EQ(vals.rank(), 2);
    EXPECT_EQ(vals.extent(0), dom_x.size());
    EXPECT_EQ(vals.extent(1), dom_y.size());

    // Fill the span with values on GPU
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            dom_x_y,
            KOKKOS_LAMBDA(IdxXY idx_x_y) {
                Idx<GridX> idx_x(idx_x_y);
                Idx<GridY> idx_y(idx_x_y);
                vals(idx_x - dom_x.front(), idx_y - dom_y.front())
                        = idx_x.uid() * nelems_y.value() + idx_y.uid();
            });

    // Get an mdspan describing the 1st and 2nd x-derivatives of the function
    auto dx_vals = dxdyField.get_mdspan(x_deriv_block);
    // Check the shape of the span
    EXPECT_EQ(dx_vals.rank(), 3);
    EXPECT_EQ(dx_vals.extent(0), x_deriv_block.size());
    EXPECT_EQ(dx_vals.extent(1), deriv_idx_range_x.size());
    EXPECT_EQ(dx_vals.extent(2), dom_y.size());

    // Collapse the index ranges into 1 iterable
    IdxRange<dX, GridY> dom_dx_y(x_deriv_block, dom_y);

    // Fill the span with values on GPU
    for (auto idx_x : deriv_idx_range_x) {
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                dom_dx_y,
                KOKKOS_LAMBDA(Idx<dX, GridY> idx_dx_y) {
                    Idx<dX> idx_dx(idx_dx_y);
                    Idx<GridY> idx_y(idx_dx_y);
                    dx_vals(idx_dx - x_deriv_block.front(),
                            deriv_idx_range_x.get_index(idx_x),
                            idx_y - dom_y.front())
                            = 100 + idx_dx.uid() * 50 + idx_x.uid() * nelems_y.value()
                              + idx_y.uid();
                });
    }

    // Get an mdspan describing the 1st and 2nd y-derivatives of the function
    auto dy_vals = dxdyField.get_mdspan(y_deriv_block);
    // Check the shape of the span
    EXPECT_EQ(dy_vals.rank(), 3);
    EXPECT_EQ(dy_vals.extent(0), y_deriv_block.size());
    EXPECT_EQ(dy_vals.extent(1), dom_x.size());
    EXPECT_EQ(dy_vals.extent(2), deriv_idx_range_y.size());

    // Collapse the index ranges into 1 iterable
    IdxRange<dY, GridX> dom_x_dy(y_deriv_block, dom_x);

    // Fill the span with values on GPU
    for (auto idx_y : deriv_idx_range_y) {
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                dom_x_dy,
                KOKKOS_LAMBDA(Idx<dY, GridX> idx_x_dy) {
                    Idx<dY> idx_dy(idx_x_dy);
                    Idx<GridX> idx_x(idx_x_dy);
                    dy_vals(idx_dy - y_deriv_block.front(),
                            idx_x - dom_x.front(),
                            deriv_idx_range_y.get_index(idx_y))
                            = 200 + idx_dy.uid() * 50 + idx_x.uid() * nelems_y.value()
                              + idx_y.uid();
                });
    }

    // Get an mdspan describing the cross-derivatives (d_xd_y, d_xd_y^2, d_x^2d_y, d_x^2d_y) of the function
    auto dx_dy_vals = dxdyField.get_mdspan(x_y_deriv_block);
    // Check the shape of the span
    EXPECT_EQ(dx_dy_vals.rank(), 4);
    EXPECT_EQ(dx_dy_vals.extent(0), x_deriv_block.size());
    EXPECT_EQ(dx_dy_vals.extent(1), y_deriv_block.size());
    EXPECT_EQ(dx_dy_vals.extent(2), deriv_idx_range_x.size());
    EXPECT_EQ(dx_dy_vals.extent(3), deriv_idx_range_y.size());

    // Collapse the index ranges into 1 iterable
    IdxRange<dX, dY> dom_dx_dy(x_deriv_block, y_deriv_block);

    // Fill the span with values on GPU
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            dom_dx_dy,
            KOKKOS_LAMBDA(Idx<dX, dY> idx_dx_dy) {
                Idx<dX> idx_dx(idx_dx_dy);
                Idx<dY> idx_dy(idx_dx_dy);
                for (auto idx_x : deriv_idx_range_x) {
                    for (auto idx_y : deriv_idx_range_y) {
                        dx_dy_vals(
                                idx_dx - x_deriv_block.front(),
                                idx_dy - y_deriv_block.front(),
                                deriv_idx_range_x.get_index(idx_x),
                                deriv_idx_range_y.get_index(idx_y))
                                = 300 + idx_dx.uid() * 100 + idx_dy.uid() * 50
                                  + idx_x.uid() * nelems_y.value() + idx_y.uid();
                    }
                }
            });

    // Copy to non-strided layout to be able to call create_mirror_view_and_copy
    FieldMem<int, IdxRangeXY> ddc_vals_alloc(dom_x_y);
    ddc::ChunkSpan ddc_vals = get_field(ddc_vals_alloc);
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            dom_x_y,
            KOKKOS_LAMBDA(IdxXY idx_x_y) { ddc_vals(idx_x_y) = dxdyField(idx_x_y); });

    // Check values
    auto vals_host = ddc::create_mirror_view_and_copy(get_const_field(ddc_vals));
    ddc::for_each(dom_x_y, [&](auto idx_x_y) {
        EXPECT_EQ(
                vals_host(idx_x_y),
                ddc::uid<GridX>(idx_x_y) * nelems_y.value() + ddc::uid<GridY>(idx_x_y));
    });

    IdxRange<GridX> dom_x_slice(dom_x.front(), IdxStep<GridX>(1));
    IdxRange<dX, GridX, GridY> dom_dx_x_y_slice(x_deriv_block, dom_x_slice, dom_y);
    // Copy to non-strided layout to be able to call create_mirror_view_and_copy
    FieldMem<int, IdxRange<dX, GridX, GridY>> ddc_dx_vals_alloc(dom_dx_x_y_slice);
    ddc::ChunkSpan ddc_dx_vals = get_field(ddc_dx_vals_alloc);
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            dom_dx_x_y_slice,
            KOKKOS_LAMBDA(Idx<dX, GridX, GridY> idx) { ddc_dx_vals(idx) = dxdyField(idx); });

    // Check x-derivative values
    auto dx_vals_host = ddc::create_mirror_view_and_copy(get_const_field(ddc_dx_vals));
    ddc::for_each(dom_dx_x_y_slice, [&](auto idx_dx_x_y) {
        Idx<dX> idx_dx(idx_dx_x_y);
        Idx<GridX> idx_x(idx_dx_x_y);
        Idx<GridY> idx_y(idx_dx_x_y);
        EXPECT_EQ(
                dx_vals_host(idx_dx_x_y),
                100 + idx_dx.uid() * 50 + idx_x.uid() * nelems_y.value() + idx_y.uid());
    });

    IdxRange<GridY> dom_y_slice(dom_y.front(), IdxStep<GridY>(1));
    IdxRange<dY, GridX, GridY> dom_dy_x_y_slice(y_deriv_block, dom_x, dom_y_slice);
    // Copy to non-strided layout to be able to call create_mirror_view_and_copy
    FieldMem<int, IdxRange<dY, GridX, GridY>> ddc_dy_vals_alloc(dom_dy_x_y_slice);
    ddc::ChunkSpan ddc_dy_vals = get_field(ddc_dy_vals_alloc);
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            dom_dy_x_y_slice,
            KOKKOS_LAMBDA(Idx<dY, GridX, GridY> idx) { ddc_dy_vals(idx) = dxdyField(idx); });

    // Check y-derivative values
    auto dy_vals_host = ddc::create_mirror_view_and_copy(get_const_field(ddc_dy_vals));
    ddc::for_each(dom_dy_x_y_slice, [&](auto idx_dy_x_y) {
        Idx<dY> idx_dy(idx_dy_x_y);
        Idx<GridX> idx_x(idx_dy_x_y);
        Idx<GridY> idx_y(idx_dy_x_y);
        EXPECT_EQ(
                dy_vals_host(idx_dy_x_y),
                200 + idx_dy.uid() * 50 + idx_x.uid() * nelems_y.value() + idx_y.uid());
    });

    IdxRange<dX, dY, GridX, GridY>
            dom_dx_dy_x_y_slice(x_deriv_block, y_deriv_block, dom_x_slice, dom_y_slice);
    // Copy to non-strided layout to be able to call create_mirror_view_and_copy
    FieldMem<int, IdxRange<dX, dY, GridX, GridY>> ddc_dx_dy_vals_alloc(dom_dx_dy_x_y_slice);
    ddc::ChunkSpan ddc_dx_dy_vals = get_field(ddc_dx_dy_vals_alloc);
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            dom_dx_dy_x_y_slice,
            KOKKOS_LAMBDA(Idx<dX, dY, GridX, GridY> idx) { ddc_dx_dy_vals(idx) = dxdyField(idx); });

    // Check cross-derivative values
    auto dx_dy_vals_host = ddc::create_mirror_view_and_copy(get_const_field(ddc_dx_dy_vals));
    ddc::for_each(dom_dx_dy_x_y_slice, [&](auto idx_dx_dy_x_y) {
        Idx<dX> idx_dx(idx_dx_dy_x_y);
        Idx<dY> idx_dy(idx_dx_dy_x_y);
        Idx<GridX> idx_x(idx_dx_dy_x_y);
        Idx<GridY> idx_y(idx_dx_dy_x_y);
        EXPECT_EQ(
                dx_dy_vals_host(idx_dx_dy_x_y),
                300 + idx_dx.uid() * 100 + idx_dy.uid() * 50 + idx_x.uid() * nelems_y.value()
                        + idx_y.uid());
    });
}

TEST(DerivFieldSpanTest, GPUElementAccess)
{
    test_DerivField_GPUElementAccess();
}

// Test if the deepcopy correctly fills a field
TEST(DerivFieldSpanTest, FieldDeepCopy)
{
    // Domains where derivatives are defined
    DiscreteSubDomain<GridX> deriv_idx_range_x(dom_x.front(), IdxStepX(2), dom_x.extents());

    // Define fields on x-y with 1 derivative in x and y
    DerivField<double, IdxRange<dX, GridX, GridY>, 1> dxdyField(dom_x_y, deriv_idx_range_x);
    DerivField<double, IdxRange<dX, GridX, GridY>, 1> dxdyField_copy(dom_x_y, deriv_idx_range_x);

    // Extract the values and derivatives
    Idx<dX> first_deriv(1);
    ddc::ChunkSpan values = dxdyField.get_values_span();
    ddc::ChunkSpan left_x_derivs = dxdyField[IdxXdX(first_deriv, deriv_idx_range_x.front())];
    ddc::ChunkSpan right_x_derivs = dxdyField[IdxXdX(first_deriv, deriv_idx_range_x.back())];

    // Set the values
    ddc::for_each(get_idx_range(values), [&](IdxXY ixy) {
        IdxX ix(ixy);
        IdxY iy(ixy);
        double x(ix - IdxX(0));
        double y(0.1 * (iy - IdxY(0)));
        values(ixy) = x * x + y * y;
    });
    // Set the derivatives
    double x_left(dom_x.front() - IdxX(0));
    ddc::for_each(get_idx_range(left_x_derivs), [&](IdxY iy) { left_x_derivs(iy) = 2 * x_left; });
    double x_right(dom_x.back() - IdxX(0));
    ddc::for_each(get_idx_range(right_x_derivs), [&](IdxY iy) {
        right_x_derivs(iy) = 2 * x_right;
    });

    // Copy the deriv field
    ddcHelper::deepcopy(dxdyField_copy, dxdyField);

    // Extract the values and derivatives from the copy
    ddc::ChunkSpan values_copy = dxdyField_copy.get_values_span();
    ddc::ChunkSpan left_x_derivs_copy = dxdyField[IdxXdX(first_deriv, deriv_idx_range_x.front())];
    ddc::ChunkSpan right_x_derivs_copy = dxdyField[IdxXdX(first_deriv, deriv_idx_range_x.back())];

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

// Test if the deepcopy correctly fills a span
TEST(DerivFieldSpanTest, SpanDeepCopy)
{
    // Domains where derivatives are defined
    DiscreteSubDomain<GridX> deriv_idx_range_x(dom_x.front(), IdxStepX(2), dom_x.extents());

    // Define fields on x-y with 1 derivative in x and y
    DerivField<double, IdxRange<dX, GridX, GridY>, 1> dxdyField_alloc(dom_x_y, deriv_idx_range_x);
    DerivField<double, IdxRange<dX, GridX, GridY>, 1>
            dxdyField_copy_alloc(dom_x_y, deriv_idx_range_x);

    // Get the spans
    DerivFieldSpan dxdyField = get_field(dxdyField_alloc);
    DerivFieldSpan dxdyField_copy = get_field(dxdyField_copy_alloc);

    // Extract the values and derivatives
    Idx<dX> first_deriv(1);
    ddc::ChunkSpan values = dxdyField.get_values_span();
    ddc::ChunkSpan left_x_derivs = dxdyField[IdxXdX(first_deriv, deriv_idx_range_x.front())];
    ddc::ChunkSpan right_x_derivs = dxdyField[IdxXdX(first_deriv, deriv_idx_range_x.back())];

    // Set the values
    ddc::for_each(get_idx_range(values), [&](IdxXY ixy) {
        IdxX ix(ixy);
        IdxY iy(ixy);
        double x(ix - IdxX(0));
        double y(0.1 * (iy - IdxY(0)));
        values(ixy) = x * x + y * y;
    });
    // Set the derivatives
    double x_left(dom_x.front() - IdxX(0));
    ddc::for_each(get_idx_range(left_x_derivs), [&](IdxY iy) { left_x_derivs(iy) = 2 * x_left; });
    double x_right(dom_x.back() - IdxX(0));
    ddc::for_each(get_idx_range(right_x_derivs), [&](IdxY iy) {
        right_x_derivs(iy) = 2 * x_right;
    });

    // Copy the deriv field
    ddcHelper::deepcopy(dxdyField_copy, dxdyField);

    // Extract the values and derivatives from the copy
    ddc::ChunkSpan values_copy = dxdyField_copy.get_values_span();
    ddc::ChunkSpan left_x_derivs_copy = dxdyField[IdxXdX(first_deriv, deriv_idx_range_x.front())];
    ddc::ChunkSpan right_x_derivs_copy = dxdyField[IdxXdX(first_deriv, deriv_idx_range_x.back())];

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
