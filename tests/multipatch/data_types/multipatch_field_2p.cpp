// SPDX-License-Identifier: MIT

/**
 * Test the MultipatchField class on two patches.
 */


#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <gtest/gtest.h>

#include "2patches_2d_non_periodic_uniform.hpp"
#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "multipatch_field.hpp"
#include "multipatch_type.hpp"
#include "patch.hpp"
#include "types.hpp"

// Namespace of the multipatch geometry where the patches are defined
using namespace non_periodic_uniform_2d_2patches;

namespace {
using DFieldMem1 = DFieldMem<Patch1::IdxRange12>;
using DFieldMem2 = DFieldMem<Patch2::IdxRange12>;

using Field1 = DFieldOnPatch<Patch1>;
using Field2 = DFieldOnPatch<Patch2>;

using ConstField1 = DConstFieldOnPatch<Patch1>;
using ConstField2 = DConstFieldOnPatch<Patch2>;

} // namespace

class MultiPatchField2Patches : public ::testing::Test
{
protected:
    // Grid sizes and intervals of the two patches
    static constexpr Patch1::IdxStep1 x1_size = Patch1::IdxStep1(5);
    static constexpr Patch1::IdxStep2 y1_size = Patch1::IdxStep2(7);

    static constexpr Patch2::IdxStep1 x2_size = Patch2::IdxStep1(9);
    static constexpr Patch2::IdxStep2 y2_size = Patch2::IdxStep2(11);

    static constexpr Patch1::Coord1 x1_min = Patch1::Coord1(0.0);
    static constexpr Patch1::Coord1 x1_max = Patch1::Coord1(1.0);

    static constexpr Patch1::Coord2 y1_min = Patch1::Coord2(0.0);
    static constexpr Patch1::Coord2 y1_max = Patch1::Coord2(1.0);

    static constexpr Patch2::Coord1 x2_min = Patch2::Coord1(0.0);
    static constexpr Patch2::Coord1 x2_max = Patch2::Coord1(1.0);

    static constexpr Patch2::Coord2 y2_min = Patch2::Coord2(0.0);
    static constexpr Patch2::Coord2 y2_max = Patch2::Coord2(1.0);


    // Index ranges of the grids on the two patches
    Patch1::IdxRange12 const index_range_xy1;
    Patch2::IdxRange12 const index_range_xy2;

    // // Memory allocations for the a field on each patch.
    DFieldMem1 field_mem1;
    DFieldMem2 field_mem2;


public:
    static void SetUpTestSuite()
    {
        // Initialising the grid on each patch

        ddc::init_discrete_space<Patch1::Grid1>(Patch1::Grid1::init(x1_min, x1_max, x1_size));
        ddc::init_discrete_space<Patch1::Grid2>(Patch1::Grid2::init(y1_min, y1_max, y1_size));

        ddc::init_discrete_space<Patch2::Grid1>(Patch2::Grid1::init(x2_min, x2_max, x2_size));
        ddc::init_discrete_space<Patch2::Grid2>(Patch2::Grid2::init(y2_min, y2_max, y2_size));
    }

    /// Fill field1 with function x^2 + y^2 and field2 with x^2 - y^2
    void initialise_non_const_fields(Field1 field1, Field2 field2)
    {
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                get_idx_range(field1),
                KOKKOS_LAMBDA(Patch1::Idx12 const idx) {
                    Patch1::Coord12 const xy1 = ddc::coordinate(idx);
                    double const x1 = ddc::get<X<1>>(xy1);
                    double const y1 = ddc::get<Y<1>>(xy1);
                    field1(idx) = x1 * x1 + y1 * y1;
                });
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                get_idx_range(field2),
                KOKKOS_LAMBDA(Patch2::Idx12 const idx) {
                    Patch2::Coord12 const xy2 = ddc::coordinate(idx);
                    double const x2 = ddc::get<X<2>>(xy2);
                    double const y2 = ddc::get<Y<2>>(xy2);
                    field2(idx) = x2 * x2 - y2 * y2;
                });
    }

    MultiPatchField2Patches()
        : index_range_xy1(
                Patch1::IdxRange1(Patch1::Idx1(0), x1_size),
                Patch1::IdxRange2(Patch1::Idx2(0), y1_size))
        , index_range_xy2(
                  Patch2::IdxRange1(Patch2::Idx1(0), x2_size),
                  Patch2::IdxRange2(Patch2::Idx2(0), y2_size))
        , field_mem1(index_range_xy1)
        , field_mem2(index_range_xy2) {};
};


/**
 * Set the field on the first patch to constant zero and the field on the 
 * second patch to constant one.
 */
TEST_F(MultiPatchField2Patches, ZeroOnFirstPatchOneOnSecond)
{
    // Arrange
    Field1 const field1 = get_field(field_mem1);
    Field2 const field2 = get_field(field_mem2);

    ddc::parallel_fill(Kokkos::DefaultExecutionSpace(), field1, 0.0);
    ddc::parallel_fill(Kokkos::DefaultExecutionSpace(), field2, 1.0);

    MultipatchField<DFieldOnPatch, Patch1, Patch2> global_field(field1, field2);

    // Act
    Field1 field1_from_multipatch = global_field.get<Patch1>();
    Field2 field2_from_multipatch = global_field.get<Patch2>();

    // Assert
    auto field1_from_multipatch_host = ddc::create_mirror_and_copy(field1_from_multipatch);
    auto field2_from_multipatch_host = ddc::create_mirror_and_copy(field2_from_multipatch);

    EXPECT_NEAR(field1_from_multipatch_host(Patch1::Idx12(4, 5)), 0.0, 1e-14);
    EXPECT_NEAR(field2_from_multipatch_host(Patch2::Idx12(7, 2)), 1.0, 1e-14);
}

/**
 * Take the scalar field x^2 + y^2 on the first patch and x^2 - y^2 on the 
 * second.
 */
TEST_F(MultiPatchField2Patches, Polynomials)
{
    // Arrange
    Field1 const field1 = get_field(field_mem1);
    Field2 const field2 = get_field(field_mem2);

    initialise_non_const_fields(field1, field2);

    MultipatchField<DFieldOnPatch, Patch1, Patch2> global_field(field1, field2);

    // Act
    Field1 field1_from_multipatch = global_field.get<Patch1>();
    Field2 field2_from_multipatch = global_field.get<Patch2>();

    // Assert
    auto field1_from_multipatch_host = ddc::create_mirror_and_copy(field1_from_multipatch);
    auto field2_from_multipatch_host = ddc::create_mirror_and_copy(field2_from_multipatch);

    std::size_t const idx_x1 = 3, idx_y1 = 5, idx_x2 = 7, idx_y2 = 2;
    double x1_length = ddc::get<Patch1::Dim1>(x1_max) - ddc::get<Patch1::Dim1>(x1_min);
    double y1_length = ddc::get<Patch1::Dim2>(y1_max) - ddc::get<Patch1::Dim2>(y1_min);
    double x2_length = ddc::get<Patch2::Dim1>(x2_max) - ddc::get<Patch2::Dim1>(x2_min);
    double y2_length = ddc::get<Patch2::Dim2>(y2_max) - ddc::get<Patch2::Dim2>(y2_min);

    EXPECT_NEAR(
            field1_from_multipatch_host(Patch1::Idx12(idx_x1, idx_y1)),
            (idx_x1 * x1_length / (x1_size.value() - 1))
                            * (idx_x1 * x1_length / (x1_size.value() - 1))
                    + (idx_y1 * y1_length / (y1_size.value() - 1))
                              * (idx_y1 * y1_length / (y1_size.value() - 1)),
            1e-14);
    EXPECT_NEAR(
            field2_from_multipatch_host(Patch2::Idx12(idx_x2, idx_y2)),
            (idx_x2 * x2_length / (x2_size.value() - 1))
                            * (idx_x2 * x2_length / (x2_size.value() - 1))
                    - (idx_y2 * y2_length / (y2_size.value() - 1))
                              * (idx_y2 * y2_length / (y2_size.value() - 1)),
            1e-14);
}

TEST_F(MultiPatchField2Patches, GetIdxRange)
{
    // Arrange
    Field1 const field1 = get_field(field_mem1);
    Field2 const field2 = get_field(field_mem2);

    initialise_non_const_fields(field1, field2);

    MultipatchField<DFieldOnPatch, Patch1, Patch2> global_field(field1, field2);

    MultipatchType<IdxRangeOnPatch, Patch1, Patch2> all_idx_ranges(get_idx_range(global_field));

    ASSERT_EQ(get_idx_range(field1), all_idx_ranges.get<Patch1>());
    ASSERT_EQ(get_idx_range(field2), all_idx_ranges.get<Patch2>());
}

TEST_F(MultiPatchField2Patches, ConstField)
{
    // Arrange
    Field1 const field1 = get_field(field_mem1);
    Field2 const field2 = get_field(field_mem2);

    initialise_non_const_fields(field1, field2);

    MultipatchField<DFieldOnPatch, Patch1, Patch2> global_field(field1, field2);

    MultipatchField<DConstFieldOnPatch, Patch1, Patch2> global_const_field(
            global_field.get_const_field());

    // Act
    ConstField1 field1_from_multipatch = global_const_field.get<Patch1>();
    ConstField2 field2_from_multipatch = global_const_field.get<Patch2>();

    // Assert
    auto field1_from_multipatch_host = ddc::create_mirror_and_copy(field1_from_multipatch);
    auto field2_from_multipatch_host = ddc::create_mirror_and_copy(field2_from_multipatch);

    std::size_t const idx_x1 = 3, idx_y1 = 5, idx_x2 = 7, idx_y2 = 2;
    double x1_length = ddc::get<Patch1::Dim1>(x1_max) - ddc::get<Patch1::Dim1>(x1_min);
    double y1_length = ddc::get<Patch1::Dim2>(y1_max) - ddc::get<Patch1::Dim2>(y1_min);
    double x2_length = ddc::get<Patch2::Dim1>(x2_max) - ddc::get<Patch2::Dim1>(x2_min);
    double y2_length = ddc::get<Patch2::Dim2>(y2_max) - ddc::get<Patch2::Dim2>(y2_min);

    EXPECT_NEAR(
            field1_from_multipatch_host(Patch1::Idx12(idx_x1, idx_y1)),
            (idx_x1 * x1_length / (x1_size.value() - 1))
                            * (idx_x1 * x1_length / (x1_size.value() - 1))
                    + (idx_y1 * y1_length / (y1_size.value() - 1))
                              * (idx_y1 * y1_length / (y1_size.value() - 1)),
            1e-14);
    EXPECT_NEAR(
            field2_from_multipatch_host(Patch2::Idx12(idx_x2, idx_y2)),
            (idx_x2 * x2_length / (x2_size.value() - 1))
                            * (idx_x2 * x2_length / (x2_size.value() - 1))
                    - (idx_y2 * y2_length / (y2_size.value() - 1))
                              * (idx_y2 * y2_length / (y2_size.value() - 1)),
            1e-14);
}
