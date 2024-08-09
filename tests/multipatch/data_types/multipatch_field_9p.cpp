// SPDX-License-Identifier: MIT

/**
 * Test the MultipatchType class by choosing three patches of a 9 patch domain
 * and defining a field on each.
 */


#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include "9patches_2d_periodic_strips_uniform.hpp"
#include "ddc_aliases.hpp"
#include "multipatch_type.hpp"

// Namespace of the multipatch geometry where the patches are defined
using namespace periodic_strips_uniform_2d_9patches;


namespace {
using DFieldMem3 = DFieldMem<Patch3::IdxRange12>;
using DFieldMem5 = DFieldMem<Patch5::IdxRange12>;
using DFieldMem7 = DFieldMem<Patch7::IdxRange12>;

template <class Patch>
using DFieldOnPatch = DField<typename Patch::IdxRange12>;

using Field3 = DFieldOnPatch<Patch3>;
using Field5 = DFieldOnPatch<Patch5>;
using Field7 = DFieldOnPatch<Patch7>;

} // namespace


TEST(MultipatchFieldMorePatches, ConstantOnEveryOne)
{
    // Arrange

    // Grid sizes of the three patches where we define fields
    static constexpr Patch3::IdxStep1 x3_size = Patch3::IdxStep1(5);
    static constexpr Patch3::IdxStep2 y3_size = Patch3::IdxStep2(7);

    static constexpr Patch5::IdxStep1 x5_size = Patch5::IdxStep1(9);
    static constexpr Patch5::IdxStep2 y5_size = Patch5::IdxStep2(11);

    static constexpr Patch7::IdxStep1 x7_size = Patch7::IdxStep1(4);
    static constexpr Patch7::IdxStep2 y7_size = Patch7::IdxStep2(8);
    Patch3::Coord1 const x3_min(0.0);
    Patch3::Coord1 const x3_max(1.0);
    Patch3::Coord2 const y3_min(0.0);
    Patch3::Coord2 const y3_max(1.0);

    Patch5::Coord1 const x5_min(0.0);
    Patch5::Coord1 const x5_max(1.0);
    Patch5::Coord2 y5_min(0.0);
    Patch5::Coord2 y5_max(1.0);

    Patch7::Coord1 const x7_min(0.0);
    Patch7::Coord1 const x7_max(1.0);
    Patch7::Coord2 y7_min(0.0);
    Patch7::Coord2 y7_max(1.0);

    ddc::init_discrete_space<Patch3::Grid1>(Patch3::Grid1::init(x3_min, x3_max, x3_size));
    ddc::init_discrete_space<Patch3::Grid2>(Patch3::Grid2::init(y3_min, y3_max, y3_size));
    ddc::init_discrete_space<Patch5::Grid1>(Patch5::Grid1::init(x5_min, x5_max, x5_size));
    ddc::init_discrete_space<Patch5::Grid2>(Patch5::Grid2::init(y5_min, y5_max, y5_size));
    ddc::init_discrete_space<Patch7::Grid1>(Patch7::Grid1::init(x7_min, x7_max, x7_size));
    ddc::init_discrete_space<Patch7::Grid2>(Patch7::Grid2::init(y7_min, y7_max, y7_size));

    Patch3::IdxRange12 const index_range_xy3(
            Patch3::IdxRange1(Patch3::Idx1(0), x3_size),
            Patch3::IdxRange2(Patch3::Idx2(0), y3_size));
    Patch5::IdxRange12 const index_range_xy5(
            Patch5::IdxRange1(Patch5::Idx1(0), x5_size),
            Patch5::IdxRange2(Patch5::Idx2(0), y5_size));
    Patch7::IdxRange12 const index_range_xy7(
            Patch7::IdxRange1(Patch7::Idx1(0), x7_size),
            Patch7::IdxRange2(Patch7::Idx2(0), y7_size));
    DFieldMem3 field_mem3(index_range_xy3);
    DFieldMem5 field_mem5(index_range_xy5);
    DFieldMem7 field_mem7(index_range_xy7);

    Field3 const field3 = get_field(field_mem3);
    Field5 const field5 = get_field(field_mem5);
    Field7 const field7 = get_field(field_mem7);

    // Set the value of field on patch k to constant k.
    ddc::parallel_fill(Kokkos::DefaultExecutionSpace(), field3, 3.0);
    ddc::parallel_fill(Kokkos::DefaultExecutionSpace(), field5, 5.0);
    ddc::parallel_fill(Kokkos::DefaultExecutionSpace(), field7, 7.0);

    MultipatchType<DFieldOnPatch, Patch3, Patch5, Patch7> global_field(field3, field5, field7);

    // Act
    Field3 field3_from_multipatch = global_field.get<Patch3>();
    Field5 field5_from_multipatch = global_field.get<Patch5>();
    Field7 field7_from_multipatch = global_field.get<Patch7>();

    // Assert
    auto field3_from_multipatch_host = ddc::create_mirror_and_copy(field3_from_multipatch);
    auto field5_from_multipatch_host = ddc::create_mirror_and_copy(field5_from_multipatch);
    auto field7_from_multipatch_host = ddc::create_mirror_and_copy(field7_from_multipatch);

    EXPECT_NEAR(field3_from_multipatch_host(Patch3::Idx12(4, 5)), 3.0, 1e-14);
    EXPECT_NEAR(field5_from_multipatch_host(Patch5::Idx12(7, 2)), 5.0, 1e-14);
    EXPECT_NEAR(field7_from_multipatch_host(Patch7::Idx12(2, 6)), 7.0, 1e-14);
}