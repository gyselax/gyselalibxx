// SPDX-License-Identifier: MIT
#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include "2patches_2d_onion_shape_uniform.hpp"
#include "ddc_aliases.hpp"
#include "multipatch_field.hpp"
#include "types.hpp"

// Namespace of the multipatch geometry where the patches are defined
using namespace onion_shape_uniform_2d_2patches;

namespace {

using VectorFieldMem1 = DVectorFieldMemOnPatch<Patch1>;
using VectorFieldMem2 = DVectorFieldMemOnPatch<Patch2>;

using Field1 = DFieldOnPatch<Patch1>;
using Field2 = DFieldOnPatch<Patch2>;

using VectorField1 = DVectorFieldOnPatch<Patch1>;
using VectorField2 = DVectorFieldOnPatch<Patch2>;

} // namespace

TEST(MultiPatchField2PatchesOnion, VectorField)
{
    Coord<R> r_min(0.0);
    Coord<R> r_mid(1.0);
    Coord<R> r_max(3.0);
    IdxStep<GridR<1>> r_size_1(5);
    IdxStep<GridR<2>> r_size_2(20);

    ddc::init_discrete_space<Patch1::Grid1>(Patch1::Grid1::init(r_min, r_mid, r_size_1));
    ddc::init_discrete_space<Patch2::Grid1>(Patch2::Grid1::init(r_mid, r_max, r_size_2));

    Coord<Theta> theta_min(0.0);
    Coord<Theta> theta_mid(M_PI);
    Coord<Theta> theta_max(2 * M_PI);
    IdxStep<GridTheta<1>> theta_size_1(10);
    IdxStep<GridTheta<2>> theta_size_2(5);

    ddc::init_discrete_space<Patch1::Grid2>(
            Patch1::Grid2::init(theta_min, theta_mid, theta_size_1));
    ddc::init_discrete_space<Patch2::Grid2>(
            Patch2::Grid2::init(theta_mid, theta_max, theta_size_2));

    IdxRange<GridR<1>> idx_range_r_1(Idx<GridR<1>> {0}, r_size_1);
    IdxRange<GridR<2>> idx_range_r_2(Idx<GridR<2>> {0}, r_size_2);

    IdxRange<GridTheta<1>> idx_range_theta_1(Idx<GridTheta<1>> {0}, theta_size_1);
    IdxRange<GridTheta<2>> idx_range_theta_2(Idx<GridTheta<2>> {0}, theta_size_2);

    IdxRange<GridR<1>, GridTheta<1>> idx_range_1(idx_range_r_1, idx_range_theta_1);
    IdxRange<GridR<2>, GridTheta<2>> idx_range_2(idx_range_r_2, idx_range_theta_2);

    // Arrange
    VectorFieldMem1 vec_field_mem1(idx_range_1);
    VectorFieldMem2 vec_field_mem2(idx_range_2);

    VectorField1 vec_field1 = get_field(vec_field_mem1);
    VectorField2 vec_field2 = get_field(vec_field_mem2);

    MultipatchField<DVectorFieldOnPatch, Patch1, Patch2> global_field(vec_field1, vec_field2);

    MultipatchField<DFieldOnPatch, Patch1, Patch2> global_field_r = ddcHelper::get<R>(global_field);
    MultipatchField<DFieldOnPatch, Patch1, Patch2> global_field_theta
            = ddcHelper::get<Theta>(global_field);

    EXPECT_EQ(
            global_field_r.template get<Patch1>().data_handle(),
            ddcHelper::get<R>(vec_field1).data_handle());
    EXPECT_EQ(
            global_field_r.template get<Patch2>().data_handle(),
            ddcHelper::get<R>(vec_field2).data_handle());
    EXPECT_EQ(
            global_field_theta.template get<Patch1>().data_handle(),
            ddcHelper::get<Theta>(vec_field1).data_handle());
    EXPECT_EQ(
            global_field_theta.template get<Patch2>().data_handle(),
            ddcHelper::get<Theta>(vec_field2).data_handle());
}
