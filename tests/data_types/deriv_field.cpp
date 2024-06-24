// SPDX-License-Identifier: MIT
#include <ddc/ddc.hpp>
#include <ddc/kernels/splines/deriv.hpp>

#include <gtest/gtest.h>

#include <ddc_helper.hpp>
#include <derivative_field.hpp>
#include <derivative_field_span.hpp>
#include <directional_tag.hpp>

namespace {

struct IDimX
{
};
using dX = ddc::Deriv<IDimX>;
using IndexX = ddc::DiscreteElement<IDimX>;
using IVectX = ddc::DiscreteVector<IDimX>;
using IDomainX = ddc::DiscreteDomain<IDimX>;


struct IDimY
{
};
using dY = ddc::Deriv<IDimY>;
using IndexY = ddc::DiscreteElement<IDimY>;
using IVectY = ddc::DiscreteVector<IDimY>;
using IDomainY = ddc::DiscreteDomain<IDimY>;

using IndexXY = ddc::DiscreteElement<IDimX, IDimY>;
using IVectXY = ddc::DiscreteVector<IDimX, IDimY>;
using IDomainXY = ddc::DiscreteDomain<IDimX, IDimY>;

using IndexXdX = ddc::DiscreteElement<dX, IDimX>;

static IndexX constexpr lbound_x(50);
static IVectX constexpr nelems_x(3);
static IDomainX constexpr dom_x(lbound_x, nelems_x);

static IndexY constexpr lbound_y(4);
static IVectY constexpr nelems_y(12);
static IDomainY constexpr dom_y(lbound_y, nelems_y);

static IndexXY constexpr lbound_x_y {lbound_x, lbound_y};
static IVectXY constexpr nelems_x_y(nelems_x, nelems_y);
static IDomainXY constexpr dom_x_y(lbound_x_y, nelems_x_y);
} // namespace

// Test the constructor for a 2D field with derivatives in 1 direction
TEST(DerivFieldTest, Constructor1Deriv)
{
    // Type for a x,y field with 1 derivative in x
    using DField = DerivField<double, ddc::DiscreteDomain<dX, IDimX, IDimY>, 1>;

    // Domain where derivatives are defined
    DiscreteSubDomain<IDimX> deriv_dom_x(dom_x.front(), IVectX(2), dom_x.extents());

    // Define the field
    DField dxField(dom_x_y, deriv_dom_x);

    // Ensure that the internal chunk has the expected type
    constexpr bool same = std::is_same_v<
            typename DField::chunk_type,
            ddc::Chunk<double, ddc::DiscreteDomain<dX, IDimX, IDimY>>>;
    EXPECT_TRUE(same);
}

// Test the constructor for a 2D field with derivatives in 2 directions
TEST(DerivFieldTest, Constructor2Deriv)
{
    // Type for a x,y field with 1 derivative in x and 1 derivative in y
    using DField = DerivField<double, ddc::DiscreteDomain<dX, dY, IDimX, IDimY>, 1>;

    // Domains where derivatives are defined
    DiscreteSubDomain<IDimX> deriv_dom_x(dom_x.front(), IVectX(2), dom_x.extents());
    DiscreteSubDomain<IDimY> deriv_dom_y(dom_y.front(), IVectY(2), dom_y.extents());

    // Define the field
    DField dxdyField(dom_x_y, deriv_dom_x, deriv_dom_y);

    // Ensure that the internal chunk has the expected type
    bool same = std::is_same_v<
            typename DField::chunk_type,
            ddc::Chunk<double, ddc::DiscreteDomain<dX, dY, IDimX, IDimY>>>;
    EXPECT_TRUE(same);
}

// Test the span constructor for a 2D field with derivatives in 1 direction
TEST(DerivFieldSpanTest, Constructor1Deriv)
{
    // Type for a x,y field span with 1 derivative in x
    using DSpan = DerivFieldSpan<double, ddc::DiscreteDomain<dX, IDimX, IDimY>>;

    // Domain where derivatives are defined
    DiscreteSubDomain<IDimX> deriv_dom_x(dom_x.front(), IVectX(2), dom_x.extents());

    // Define a field on x-y with 1 derivative in x and y
    DerivField<double, ddc::DiscreteDomain<dX, IDimX, IDimY>, 1>
            dxField_alloc(dom_x_y, deriv_dom_x);
    // Define the span
    DSpan dxField(dxField_alloc);

    // Ensure that the internal chunk span has the expected type
    bool same = std::is_same_v<
            DSpan::chunk_type,
            ddc::ChunkSpan<double, ddc::DiscreteDomain<dX, IDimX, IDimY>>>;
    EXPECT_TRUE(same);
}

// Test the span constructor for a 2D field with derivatives in 2 directions
TEST(DerivFieldSpanTest, Constructor2Deriv)
{
    // Type for a x,y field span with 1 derivative in x and 1 derivative in y
    using DSpan = DerivFieldSpan<double, ddc::DiscreteDomain<dX, dY, IDimX, IDimY>>;

    // Domains where derivatives are defined
    DiscreteSubDomain<IDimX> deriv_dom_x(dom_x.front(), IVectX(2), dom_x.extents());
    DiscreteSubDomain<IDimY> deriv_dom_y(dom_y.front(), IVectY(2), dom_y.extents());

    // Define a field on x-y with 1 derivative in x and y
    DerivField<double, ddc::DiscreteDomain<dX, dY, IDimX, IDimY>, 1>
            dxdyField_alloc(dom_x_y, deriv_dom_x, deriv_dom_y);
    // Define the span
    DSpan dxdyField(dxdyField_alloc);

    // Ensure that the internal chunk span has the expected type
    bool same = std::is_same_v<
            DSpan::chunk_type,
            ddc::ChunkSpan<double, ddc::DiscreteDomain<dX, dY, IDimX, IDimY>>>;
    EXPECT_TRUE(same);
}

// Test if the values of the function can be accessed via the get_values_span function
TEST(DerivFieldSpanTest, SpanValueAccess)
{
    // Domains where derivatives are defined
    DiscreteSubDomain<IDimX> deriv_dom_x(dom_x.front(), IVectX(2), dom_x.extents());
    DiscreteSubDomain<IDimY> deriv_dom_y(dom_y.front(), IVectY(2), dom_y.extents());

    // Define a field on x-y with 1 derivative in x and y
    DerivField<double, ddc::DiscreteDomain<dX, dY, IDimX, IDimY>, 1>
            dxdyField_alloc(dom_x_y, deriv_dom_x, deriv_dom_y);
    // Define the span
    DerivFieldSpan<double, ddc::DiscreteDomain<dX, dY, IDimX, IDimY>> dxdyField(dxdyField_alloc);

    // Check that the domain of the values matches the expected domain
    EXPECT_EQ(dom_x_y, dxdyField.get_values_span().domain());
}

// Test if the values of the function can be accessed via the constant get_values_span function
TEST(DerivFieldSpanTest, ViewValueAccess)
{
    // Domains where derivatives are defined
    DiscreteSubDomain<IDimX> deriv_dom_x(dom_x.front(), IVectX(2), dom_x.extents());
    DiscreteSubDomain<IDimY> deriv_dom_y(dom_y.front(), IVectY(2), dom_y.extents());

    // Define a field on x-y with 1 derivative in x and y
    DerivField<double, ddc::DiscreteDomain<dX, dY, IDimX, IDimY>, 1>
            dxdyField_alloc(dom_x_y, deriv_dom_x, deriv_dom_y);
    // Define the span
    DerivFieldView<const double, ddc::DiscreteDomain<dX, dY, IDimX, IDimY>> dxdyField(
            dxdyField_alloc);

    // Check that the domain of the values matches the expected domain
    EXPECT_EQ(dom_x_y, dxdyField.get_values_span().domain());
}


// Test if the derivatives of the function can be accessed via the slice function
TEST(DerivFieldTest, derivValueAccess)
{
    // Domains where derivatives are defined
    DiscreteSubDomain<IDimX> deriv_dom_x(dom_x.front(), IVectX(2), dom_x.extents());
    DiscreteSubDomain<IDimY> deriv_dom_y(dom_y.front(), IVectY(2), dom_y.extents());

    // Define a field on x-y with 2 derivatives in x and y
    DerivField<double, ddc::DiscreteDomain<dX, dY, IDimX, IDimY>, 2>
            dxdyField(dom_x_y, deriv_dom_x, deriv_dom_y);

    // The derivatives to be retrieved from the field
    ddc::DiscreteDomain<dX> x_deriv_block(ddc::DiscreteElement<dX>(1), ddc::DiscreteVector<dX>(2));
    ddc::DiscreteDomain<dY> y_deriv_block(ddc::DiscreteElement<dY>(1), ddc::DiscreteVector<dY>(2));

    // The subdomains to be retrieved from the field
    ddc::DiscreteDomain<IDimX> deriv_x_subdom(dom_x.front(), ddc::DiscreteVector<IDimX>(1));
    ddc::DiscreteDomain<IDimY> deriv_y_subdom(dom_y.front(), ddc::DiscreteVector<IDimY>(1));

    // The expected domain for the x-derivatives identified in x_deriv_block at the positions identified by deriv_x_subdom
    ddc::DiscreteDomain<dX, IDimX, IDimY>
            dx_domain(x_deriv_block, deriv_x_subdom, ddc::select<IDimY>(dom_x_y));
    // The expected domain for the y-derivatives identified in y_deriv_block at the positions identified by deriv_y_subdom
    ddc::DiscreteDomain<dY, IDimX, IDimY>
            dy_domain(y_deriv_block, ddc::select<IDimX>(dom_x_y), deriv_y_subdom);
    // The expected domain for the x-y-derivatives identified in x_deriv_block and y_deriv_block at the positions identified
    // by deriv_x_subdom and deriv_y_subdom
    ddc::DiscreteDomain<dX, dY, IDimX, IDimY>
            dx_dy_domain(x_deriv_block, y_deriv_block, deriv_x_subdom, deriv_y_subdom);

    // Build the multi-D domains that will slice the field
    ddc::DiscreteDomain<dX, IDimX> slice_idx_dx(x_deriv_block, deriv_x_subdom);
    ddc::DiscreteDomain<dY, IDimY> slice_idx_dy(y_deriv_block, deriv_y_subdom);
    ddc::DiscreteDomain<dX, dY, IDimX, IDimY>
            slice_idx_dx_dy(x_deriv_block, y_deriv_block, deriv_x_subdom, deriv_y_subdom);

    // Collect the domain of the sliced fields
    ddc::DiscreteDomain<dX, IDimX, IDimY> slice_dx_domain = dxdyField[slice_idx_dx].domain();
    ddc::DiscreteDomain<dY, IDimX, IDimY> slice_dy_domain = dxdyField[slice_idx_dy].domain();
    ddc::DiscreteDomain<dX, dY, IDimX, IDimY> slice_dx_dy_domain
            = dxdyField[slice_idx_dx_dy].domain();

    // Check that the domains are as expected
    EXPECT_EQ(dx_domain, slice_dx_domain);
    EXPECT_EQ(dy_domain, slice_dy_domain);
    EXPECT_EQ(dx_dy_domain, slice_dx_dy_domain);
}

// Test if the derivatives of the function can be accessed from a span via the slice function
TEST(DerivFieldSpanTest, derivSpanValueAccess)
{
    // Domains where derivatives are defined
    DiscreteSubDomain<IDimX> deriv_dom_x(dom_x.front(), IVectX(2), dom_x.extents());
    DiscreteSubDomain<IDimY> deriv_dom_y(dom_y.front(), IVectY(2), dom_y.extents());

    // Define a field on x-y with 2 derivatives in x and y
    DerivField<double, ddc::DiscreteDomain<dX, dY, IDimX, IDimY>, 2>
            dxdyField_alloc(dom_x_y, deriv_dom_x, deriv_dom_y);
    // Define the span
    DerivFieldSpan<double, ddc::DiscreteDomain<dX, dY, IDimX, IDimY>> dxdyField(dxdyField_alloc);

    // The derivatives to be retrieved from the field
    ddc::DiscreteDomain<dX> x_deriv_block(ddc::DiscreteElement<dX>(1), ddc::DiscreteVector<dX>(2));
    ddc::DiscreteDomain<dY> y_deriv_block(ddc::DiscreteElement<dY>(1), ddc::DiscreteVector<dY>(2));

    // The subdomains to be retrieved from the field
    ddc::DiscreteDomain<IDimX> deriv_x_subdom(dom_x.front(), ddc::DiscreteVector<IDimX>(1));
    ddc::DiscreteDomain<IDimY> deriv_y_subdom(dom_y.front(), ddc::DiscreteVector<IDimY>(1));

    // The expected domain for the x-derivatives identified in x_deriv_block at the positions identified by deriv_x_subdom
    ddc::DiscreteDomain<dX, IDimX, IDimY>
            dx_domain(x_deriv_block, deriv_x_subdom, ddc::select<IDimY>(dom_x_y));
    // The expected domain for the y-derivatives identified in y_deriv_block at the positions identified by deriv_y_subdom
    ddc::DiscreteDomain<dY, IDimX, IDimY>
            dy_domain(y_deriv_block, ddc::select<IDimX>(dom_x_y), deriv_y_subdom);
    // The expected domain for the x-y-derivatives identified in x_deriv_block and y_deriv_block at the positions identified
    // by deriv_x_subdom and deriv_y_subdom
    ddc::DiscreteDomain<dX, dY, IDimX, IDimY>
            dx_dy_domain(x_deriv_block, y_deriv_block, deriv_x_subdom, deriv_y_subdom);

    // Build the multi-D domains that will slice the field
    ddc::DiscreteDomain<dX, IDimX> slice_idx_dx(x_deriv_block, deriv_x_subdom);
    ddc::DiscreteDomain<dY, IDimY> slice_idx_dy(y_deriv_block, deriv_y_subdom);
    ddc::DiscreteDomain<dX, dY, IDimX, IDimY>
            slice_idx_dx_dy(x_deriv_block, y_deriv_block, deriv_x_subdom, deriv_y_subdom);

    // Collect the domain of the sliced fields
    ddc::DiscreteDomain<dX, IDimX, IDimY> slice_dx_domain = dxdyField[slice_idx_dx].domain();
    ddc::DiscreteDomain<dY, IDimX, IDimY> slice_dy_domain = dxdyField[slice_idx_dy].domain();
    ddc::DiscreteDomain<dX, dY, IDimX, IDimY> slice_dx_dy_domain
            = dxdyField[slice_idx_dx_dy].domain();

    // Check that the domains are as expected
    EXPECT_EQ(dx_domain, slice_dx_domain);
    EXPECT_EQ(dy_domain, slice_dy_domain);
    EXPECT_EQ(dx_dy_domain, slice_dx_dy_domain);
}

// Test the element-wise operator
TEST(DerivFieldSpanTest, ElementAccess)
{
    // Domains where derivatives are defined
    DiscreteSubDomain<IDimX> deriv_dom_x(dom_x.front(), IVectX(2), dom_x.extents());
    DiscreteSubDomain<IDimY> deriv_dom_y(dom_y.front(), IVectY(2), dom_y.extents());

    // Define a field on x-y with 3 derivatives in x and y
    DerivField<int, ddc::DiscreteDomain<dX, dY, IDimX, IDimY>, 3>
            dxdyField_alloc(dom_x_y, deriv_dom_x, deriv_dom_y);
    // Define the span
    DerivFieldSpan<int, ddc::DiscreteDomain<dX, dY, IDimX, IDimY>> dxdyField(dxdyField_alloc);

    // A subset  of the x-derivatives to be retrieved with get_mdspan
    ddc::DiscreteDomain<dX> x_deriv_block(ddc::DiscreteElement<dX>(1), ddc::DiscreteVector<dX>(2));
    // A subset  of the y-derivatives to be retrieved with get_mdspan
    ddc::DiscreteDomain<dY> y_deriv_block(ddc::DiscreteElement<dY>(1), ddc::DiscreteVector<dY>(2));
    // A subset  of the x and y derivatives to be retrieved with get_mdspan
    ddc::DiscreteDomain<dX, dY> x_y_deriv_block(x_deriv_block, y_deriv_block);

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
    EXPECT_EQ(dx_vals.extent(1), deriv_dom_x.size());
    EXPECT_EQ(dx_vals.extent(2), dom_y.size());

    // Fill the span with values
    for (auto idx_x : deriv_dom_x) {
        ddc::for_each(dom_y, [&](auto idx_y) {
            ddc::for_each(x_deriv_block, [&](auto idx_dx) {
                dx_vals(idx_dx - x_deriv_block.front(),
                        deriv_dom_x.get_index(idx_x),
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
    EXPECT_EQ(dy_vals.extent(2), deriv_dom_y.size());

    // Fill the span with values
    for (auto idx_y : deriv_dom_y) {
        ddc::for_each(dom_x, [&](auto idx_x) {
            ddc::for_each(y_deriv_block, [&](auto idx_dy) {
                dy_vals(idx_dy - y_deriv_block.front(),
                        idx_x - dom_x.front(),
                        deriv_dom_y.get_index(idx_y))
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
    EXPECT_EQ(dx_dy_vals.extent(2), deriv_dom_x.size());
    EXPECT_EQ(dx_dy_vals.extent(3), deriv_dom_y.size());

    // Fill the span with values
    ddc::for_each(x_deriv_block, [&](auto idx_dx) {
        ddc::for_each(y_deriv_block, [&](auto idx_dy) {
            for (auto idx_x : deriv_dom_x) {
                for (auto idx_y : deriv_dom_y) {
                    dx_dy_vals(
                            idx_dx - x_deriv_block.front(),
                            idx_dy - y_deriv_block.front(),
                            deriv_dom_x.get_index(idx_x),
                            deriv_dom_y.get_index(idx_y))
                            = 300 + idx_dx.uid() * 100 + idx_dy.uid() * 50
                              + idx_x.uid() * nelems_y.value() + idx_y.uid();
                }
            }
        });
    });

    // Check that the value of the function at the index val_element is correct
    IndexXY val_element(lbound_x + 2, lbound_y + 4);
    EXPECT_EQ(
            dxdyField(val_element),
            ddc::uid<IDimX>(val_element) * nelems_y.value() + ddc::uid<IDimY>(val_element));

    // Check that the value of the 0-th derivatives of the function at the index val_element is correct
    ddc::DiscreteElement<dX, dY, IDimX, IDimY> exact_val_element(
            ddc::DiscreteElement<dX>(0),
            ddc::DiscreteElement<dY>(0),
            lbound_x + 2,
            lbound_y + 4);
    EXPECT_EQ(
            dxdyField(exact_val_element),
            ddc::uid<IDimX>(val_element) * nelems_y.value() + ddc::uid<IDimY>(val_element));

    // Check that the value of the 1st x-derivative of the function at the index (lbound_x, lbound_y) is correct
    ddc::DiscreteElement<dX, IDimX, IDimY>
            dx_element(ddc::DiscreteElement<dX>(1), lbound_x, lbound_y);
    EXPECT_EQ(
            dxdyField(dx_element),
            100 + ddc::uid<dX>(dx_element) * 50 + ddc::uid<IDimX>(dx_element) * nelems_y.value()
                    + ddc::uid<IDimY>(dx_element));

    // Check that the value of the 1st y-derivative of the function at the index (lbound_x, lbound_y) is correct
    ddc::DiscreteElement<dY, IDimX, IDimY>
            dy_element(ddc::DiscreteElement<dY>(1), lbound_x, lbound_y);
    EXPECT_EQ(
            dxdyField(dy_element),
            200 + ddc::uid<dY>(dy_element) * 50 + ddc::uid<IDimX>(dy_element) * nelems_y.value()
                    + ddc::uid<IDimY>(dy_element));

    // Check that the value of the cross derivative d_xd_y of the function at the index (lbound_x, lbound_y) is correct
    ddc::DiscreteElement<dX, dY, IDimX, IDimY> dx_dy_element(
            ddc::DiscreteElement<dX>(1),
            ddc::DiscreteElement<dY>(1),
            lbound_x,
            lbound_y);
    EXPECT_EQ(
            dxdyField(dx_dy_element),
            300 + ddc::uid<dX>(dx_dy_element) * 100 + ddc::uid<dY>(dx_dy_element) * 50
                    + ddc::uid<IDimX>(dx_dy_element) * nelems_y.value()
                    + ddc::uid<IDimY>(dx_dy_element));
}

// Test the element-wise operator on GPU
void test_DerivField_GPUElementAccess()
{
    // Domains where derivatives are defined
    DiscreteSubDomain<IDimX> deriv_dom_x(dom_x.front(), IVectX(2), dom_x.extents());
    DiscreteSubDomain<IDimY> deriv_dom_y(dom_y.front(), IVectY(2), dom_y.extents());

    // Define a field on x-y with 3 derivatives in x and y on GPU
    device_t<DerivField<int, ddc::DiscreteDomain<dX, dY, IDimX, IDimY>, 3>>
            dxdyField_alloc(dom_x_y, deriv_dom_x, deriv_dom_y);
    // Define the span on GPU
    device_t<DerivFieldSpan<int, ddc::DiscreteDomain<dX, dY, IDimX, IDimY>>> dxdyField(
            dxdyField_alloc);

    // A subset  of the x-derivatives to be retrieved with get_mdspan
    ddc::DiscreteDomain<dX> x_deriv_block(ddc::DiscreteElement<dX>(1), ddc::DiscreteVector<dX>(2));
    ddc::DiscreteDomain<dY> y_deriv_block(ddc::DiscreteElement<dY>(1), ddc::DiscreteVector<dY>(2));
    ddc::DiscreteDomain<dX, dY> x_y_deriv_block(x_deriv_block, y_deriv_block);

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
            KOKKOS_LAMBDA(IndexXY idx_x_y) {
                ddc::DiscreteElement<IDimX> idx_x(idx_x_y);
                ddc::DiscreteElement<IDimY> idx_y(idx_x_y);
                vals(idx_x - dom_x.front(), idx_y - dom_y.front())
                        = idx_x.uid() * nelems_y.value() + idx_y.uid();
            });

    // Get an mdspan describing the 1st and 2nd x-derivatives of the function
    auto dx_vals = dxdyField.get_mdspan(x_deriv_block);
    // Check the shape of the span
    EXPECT_EQ(dx_vals.rank(), 3);
    EXPECT_EQ(dx_vals.extent(0), x_deriv_block.size());
    EXPECT_EQ(dx_vals.extent(1), deriv_dom_x.size());
    EXPECT_EQ(dx_vals.extent(2), dom_y.size());

    // Collapse the domains into 1 iterable
    ddc::DiscreteDomain<dX, IDimY> dom_dx_y(x_deriv_block, dom_y);

    // Fill the span with values on GPU
    for (auto idx_x : deriv_dom_x) {
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                dom_dx_y,
                KOKKOS_LAMBDA(ddc::DiscreteElement<dX, IDimY> idx_dx_y) {
                    ddc::DiscreteElement<dX> idx_dx(idx_dx_y);
                    ddc::DiscreteElement<IDimY> idx_y(idx_dx_y);
                    dx_vals(idx_dx - x_deriv_block.front(),
                            deriv_dom_x.get_index(idx_x),
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
    EXPECT_EQ(dy_vals.extent(2), deriv_dom_y.size());

    // Collapse the domains into 1 iterable
    ddc::DiscreteDomain<dY, IDimX> dom_x_dy(y_deriv_block, dom_x);

    // Fill the span with values on GPU
    for (auto idx_y : deriv_dom_y) {
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                dom_x_dy,
                KOKKOS_LAMBDA(ddc::DiscreteElement<dY, IDimX> idx_x_dy) {
                    ddc::DiscreteElement<dY> idx_dy(idx_x_dy);
                    ddc::DiscreteElement<IDimX> idx_x(idx_x_dy);
                    dy_vals(idx_dy - y_deriv_block.front(),
                            idx_x - dom_x.front(),
                            deriv_dom_y.get_index(idx_y))
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
    EXPECT_EQ(dx_dy_vals.extent(2), deriv_dom_x.size());
    EXPECT_EQ(dx_dy_vals.extent(3), deriv_dom_y.size());

    // Collapse the domains into 1 iterable
    ddc::DiscreteDomain<dX, dY> dom_dx_dy(x_deriv_block, y_deriv_block);

    // Fill the span with values on GPU
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            dom_dx_dy,
            KOKKOS_LAMBDA(ddc::DiscreteElement<dX, dY> idx_dx_dy) {
                ddc::DiscreteElement<dX> idx_dx(idx_dx_dy);
                ddc::DiscreteElement<dY> idx_dy(idx_dx_dy);
                for (auto idx_x : deriv_dom_x) {
                    for (auto idx_y : deriv_dom_y) {
                        dx_dy_vals(
                                idx_dx - x_deriv_block.front(),
                                idx_dy - y_deriv_block.front(),
                                deriv_dom_x.get_index(idx_x),
                                deriv_dom_y.get_index(idx_y))
                                = 300 + idx_dx.uid() * 100 + idx_dy.uid() * 50
                                  + idx_x.uid() * nelems_y.value() + idx_y.uid();
                    }
                }
            });

    // Copy to non-strided layout to be able to call create_mirror_view_and_copy
    device_t<ddc::Chunk<int, IDomainXY>> ddc_vals_alloc(dom_x_y);
    ddc::ChunkSpan ddc_vals = ddc_vals_alloc.span_view();
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            dom_x_y,
            KOKKOS_LAMBDA(IndexXY idx_x_y) { ddc_vals(idx_x_y) = dxdyField(idx_x_y); });

    // Check values
    auto vals_host = ddc::create_mirror_view_and_copy(ddc_vals.span_cview());
    ddc::for_each(dom_x_y, [&](auto idx_x_y) {
        EXPECT_EQ(
                vals_host(idx_x_y),
                ddc::uid<IDimX>(idx_x_y) * nelems_y.value() + ddc::uid<IDimY>(idx_x_y));
    });

    ddc::DiscreteDomain<IDimX> dom_x_slice(dom_x.front(), ddc::DiscreteVector<IDimX>(1));
    ddc::DiscreteDomain<dX, IDimX, IDimY> dom_dx_x_y_slice(x_deriv_block, dom_x_slice, dom_y);
    // Copy to non-strided layout to be able to call create_mirror_view_and_copy
    device_t<ddc::Chunk<int, ddc::DiscreteDomain<dX, IDimX, IDimY>>> ddc_dx_vals_alloc(
            dom_dx_x_y_slice);
    ddc::ChunkSpan ddc_dx_vals = ddc_dx_vals_alloc.span_view();
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            dom_dx_x_y_slice,
            KOKKOS_LAMBDA(ddc::DiscreteElement<dX, IDimX, IDimY> idx) {
                ddc_dx_vals(idx) = dxdyField(idx);
            });

    // Check x-derivative values
    auto dx_vals_host = ddc::create_mirror_view_and_copy(ddc_dx_vals.span_cview());
    ddc::for_each(dom_dx_x_y_slice, [&](auto idx_dx_x_y) {
        ddc::DiscreteElement<dX> idx_dx(idx_dx_x_y);
        ddc::DiscreteElement<IDimX> idx_x(idx_dx_x_y);
        ddc::DiscreteElement<IDimY> idx_y(idx_dx_x_y);
        EXPECT_EQ(
                dx_vals_host(idx_dx_x_y),
                100 + idx_dx.uid() * 50 + idx_x.uid() * nelems_y.value() + idx_y.uid());
    });

    ddc::DiscreteDomain<IDimY> dom_y_slice(dom_y.front(), ddc::DiscreteVector<IDimY>(1));
    ddc::DiscreteDomain<dY, IDimX, IDimY> dom_dy_x_y_slice(y_deriv_block, dom_x, dom_y_slice);
    // Copy to non-strided layout to be able to call create_mirror_view_and_copy
    device_t<ddc::Chunk<int, ddc::DiscreteDomain<dY, IDimX, IDimY>>> ddc_dy_vals_alloc(
            dom_dy_x_y_slice);
    ddc::ChunkSpan ddc_dy_vals = ddc_dy_vals_alloc.span_view();
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            dom_dy_x_y_slice,
            KOKKOS_LAMBDA(ddc::DiscreteElement<dY, IDimX, IDimY> idx) {
                ddc_dy_vals(idx) = dxdyField(idx);
            });

    // Check y-derivative values
    auto dy_vals_host = ddc::create_mirror_view_and_copy(ddc_dy_vals.span_cview());
    ddc::for_each(dom_dy_x_y_slice, [&](auto idx_dy_x_y) {
        ddc::DiscreteElement<dY> idx_dy(idx_dy_x_y);
        ddc::DiscreteElement<IDimX> idx_x(idx_dy_x_y);
        ddc::DiscreteElement<IDimY> idx_y(idx_dy_x_y);
        EXPECT_EQ(
                dy_vals_host(idx_dy_x_y),
                200 + idx_dy.uid() * 50 + idx_x.uid() * nelems_y.value() + idx_y.uid());
    });

    ddc::DiscreteDomain<dX, dY, IDimX, IDimY>
            dom_dx_dy_x_y_slice(x_deriv_block, y_deriv_block, dom_x_slice, dom_y_slice);
    // Copy to non-strided layout to be able to call create_mirror_view_and_copy
    device_t<ddc::Chunk<int, ddc::DiscreteDomain<dX, dY, IDimX, IDimY>>> ddc_dx_dy_vals_alloc(
            dom_dx_dy_x_y_slice);
    ddc::ChunkSpan ddc_dx_dy_vals = ddc_dx_dy_vals_alloc.span_view();
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            dom_dx_dy_x_y_slice,
            KOKKOS_LAMBDA(ddc::DiscreteElement<dX, dY, IDimX, IDimY> idx) {
                ddc_dx_dy_vals(idx) = dxdyField(idx);
            });

    // Check cross-derivative values
    auto dx_dy_vals_host = ddc::create_mirror_view_and_copy(ddc_dx_dy_vals.span_cview());
    ddc::for_each(dom_dx_dy_x_y_slice, [&](auto idx_dx_dy_x_y) {
        ddc::DiscreteElement<dX> idx_dx(idx_dx_dy_x_y);
        ddc::DiscreteElement<dY> idx_dy(idx_dx_dy_x_y);
        ddc::DiscreteElement<IDimX> idx_x(idx_dx_dy_x_y);
        ddc::DiscreteElement<IDimY> idx_y(idx_dx_dy_x_y);
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
    DiscreteSubDomain<IDimX> deriv_dom_x(dom_x.front(), IVectX(2), dom_x.extents());

    // Define fields on x-y with 1 derivative in x and y
    DerivField<double, ddc::DiscreteDomain<dX, IDimX, IDimY>, 1> dxdyField(dom_x_y, deriv_dom_x);
    DerivField<double, ddc::DiscreteDomain<dX, IDimX, IDimY>, 1>
            dxdyField_copy(dom_x_y, deriv_dom_x);

    // Extract the values and derivatives
    ddc::DiscreteElement<dX> first_deriv(1);
    ddc::ChunkSpan values = dxdyField.get_values_span();
    ddc::ChunkSpan left_x_derivs = dxdyField[IndexXdX(first_deriv, deriv_dom_x.front())];
    ddc::ChunkSpan right_x_derivs = dxdyField[IndexXdX(first_deriv, deriv_dom_x.back())];

    // Set the values
    ddc::for_each(values.domain(), [&](IndexXY ixy) {
        IndexX ix(ixy);
        IndexY iy(ixy);
        double x(ix - IndexX(0));
        double y(0.1 * (iy - IndexY(0)));
        values(ixy) = x * x + y * y;
    });
    // Set the derivatives
    double x_left(dom_x.front() - IndexX(0));
    ddc::for_each(left_x_derivs.domain(), [&](IndexY iy) { left_x_derivs(iy) = 2 * x_left; });
    double x_right(dom_x.back() - IndexX(0));
    ddc::for_each(right_x_derivs.domain(), [&](IndexY iy) { right_x_derivs(iy) = 2 * x_right; });

    // Copy the deriv field
    ddcHelper::deepcopy(dxdyField_copy, dxdyField);

    // Extract the values and derivatives from the copy
    ddc::ChunkSpan values_copy = dxdyField_copy.get_values_span();
    ddc::ChunkSpan left_x_derivs_copy = dxdyField[IndexXdX(first_deriv, deriv_dom_x.front())];
    ddc::ChunkSpan right_x_derivs_copy = dxdyField[IndexXdX(first_deriv, deriv_dom_x.back())];

    // Check the values
    ddc::for_each(values.domain(), [&](IndexXY ixy) {
        IndexX ix(ixy);
        IndexY iy(ixy);
        double x(ix - IndexX(0));
        double y(0.1 * (iy - IndexY(0)));
        EXPECT_EQ(values_copy(ixy), x * x + y * y);
    });
    // Check the derivatives
    ddc::for_each(left_x_derivs.domain(), [&](IndexY iy) {
        EXPECT_EQ(left_x_derivs(iy), 2 * x_left);
    });
    ddc::for_each(right_x_derivs.domain(), [&](IndexY iy) {
        EXPECT_EQ(right_x_derivs(iy), 2 * x_right);
    });
}

// Test if the deepcopy correctly fills a span
TEST(DerivFieldSpanTest, SpanDeepCopy)
{
    // Domains where derivatives are defined
    DiscreteSubDomain<IDimX> deriv_dom_x(dom_x.front(), IVectX(2), dom_x.extents());

    // Define fields on x-y with 1 derivative in x and y
    DerivField<double, ddc::DiscreteDomain<dX, IDimX, IDimY>, 1>
            dxdyField_alloc(dom_x_y, deriv_dom_x);
    DerivField<double, ddc::DiscreteDomain<dX, IDimX, IDimY>, 1>
            dxdyField_copy_alloc(dom_x_y, deriv_dom_x);

    // Get the spans
    DerivFieldSpan dxdyField = dxdyField_alloc.span_view();
    DerivFieldSpan dxdyField_copy = dxdyField_copy_alloc.span_view();

    // Extract the values and derivatives
    ddc::DiscreteElement<dX> first_deriv(1);
    ddc::ChunkSpan values = dxdyField.get_values_span();
    ddc::ChunkSpan left_x_derivs = dxdyField[IndexXdX(first_deriv, deriv_dom_x.front())];
    ddc::ChunkSpan right_x_derivs = dxdyField[IndexXdX(first_deriv, deriv_dom_x.back())];

    // Set the values
    ddc::for_each(values.domain(), [&](IndexXY ixy) {
        IndexX ix(ixy);
        IndexY iy(ixy);
        double x(ix - IndexX(0));
        double y(0.1 * (iy - IndexY(0)));
        values(ixy) = x * x + y * y;
    });
    // Set the derivatives
    double x_left(dom_x.front() - IndexX(0));
    ddc::for_each(left_x_derivs.domain(), [&](IndexY iy) { left_x_derivs(iy) = 2 * x_left; });
    double x_right(dom_x.back() - IndexX(0));
    ddc::for_each(right_x_derivs.domain(), [&](IndexY iy) { right_x_derivs(iy) = 2 * x_right; });

    // Copy the deriv field
    ddcHelper::deepcopy(dxdyField_copy, dxdyField);

    // Extract the values and derivatives from the copy
    ddc::ChunkSpan values_copy = dxdyField_copy.get_values_span();
    ddc::ChunkSpan left_x_derivs_copy = dxdyField[IndexXdX(first_deriv, deriv_dom_x.front())];
    ddc::ChunkSpan right_x_derivs_copy = dxdyField[IndexXdX(first_deriv, deriv_dom_x.back())];

    // Check the values
    ddc::for_each(values.domain(), [&](IndexXY ixy) {
        IndexX ix(ixy);
        IndexY iy(ixy);
        double x(ix - IndexX(0));
        double y(0.1 * (iy - IndexY(0)));
        EXPECT_EQ(values(ixy), x * x + y * y);
    });
    // Check the derivatives
    ddc::for_each(left_x_derivs.domain(), [&](IndexY iy) {
        EXPECT_EQ(left_x_derivs(iy), 2 * x_left);
    });
    ddc::for_each(right_x_derivs.domain(), [&](IndexY iy) {
        EXPECT_EQ(right_x_derivs(iy), 2 * x_right);
    });
}
