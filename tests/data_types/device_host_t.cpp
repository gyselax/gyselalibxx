// SPDX-License-Identifier: MIT
#include <typeinfo>

#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
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


struct GridX
{
};
using IdxX = Idx<GridX>;
using IdxSteptX = IdxStep<GridX>;
using IdxRangeX = IdxRange<GridX>;


template <class Datatype>
using FieldMemX = FieldMem<Datatype, IdxRangeX>;
template <class Datatype>
using FieldX = Field<Datatype, IdxRangeX>;


template <class Datatype>
using VectorFieldMemX = VectorFieldMem<Datatype, IdxRangeX, Direction>;
template <class Datatype>
using VectorFieldX = VectorField<Datatype, IdxRangeX, Direction>;


static IdxX constexpr lbound_x(50);
static IdxSteptX constexpr nelems_x(3);
static IdxRangeX constexpr idx_range_x(lbound_x, nelems_x);

} // end namespace


// Tests on DefaultExecutionSpace ----------------------------------------------------------------
TEST(MemorySpace, FieldMemOnDeviceT)
{
    FieldMemX<double> field_test_alloc(idx_range_x);
    DFieldMem<IdxRangeX, ddc::KokkosAllocator<double, Kokkos::DefaultExecutionSpace::memory_space>>
            field_alloc(idx_range_x);
    EXPECT_TRUE(typeid(field_test_alloc) == typeid(field_alloc));
}


TEST(MemorySpace, FieldOnDeviceT)
{
    FieldMemX<double> field_test_alloc(idx_range_x);
    FieldX<double> field_test = get_field(field_test_alloc);
    DFieldMem<IdxRangeX, ddc::KokkosAllocator<double, Kokkos::DefaultExecutionSpace::memory_space>>
            field_alloc(idx_range_x);
    DField<IdxRangeX> field = get_field(field_alloc);
    EXPECT_TRUE(typeid(field_test) == typeid(field));
}


TEST(MemorySpace, VectorFieldMemOnDeviceT)
{
    device_t<VectorFieldMemX<double>> vector_field_test_alloc(idx_range_x);
    VectorFieldMem<
            double,
            IdxRangeX,
            Direction,
            ddc::KokkosAllocator<double, Kokkos::DefaultExecutionSpace::memory_space>>
            vector_field_alloc(idx_range_x);
    EXPECT_TRUE(typeid(vector_field_test_alloc) == typeid(vector_field_alloc));
}


TEST(MemorySpace, VectorFieldOnDeviceT)
{
    device_t<VectorFieldMemX<double>> vector_field_test_alloc(idx_range_x);
    device_t<VectorFieldX<double>> vector_field_test = get_field(vector_field_test_alloc);
    VectorFieldMem<
            double,
            IdxRangeX,
            Direction,
            ddc::KokkosAllocator<double, Kokkos::DefaultExecutionSpace::memory_space>>
            vector_field_alloc(idx_range_x);
    VectorField vector_field = get_field(vector_field_alloc);
    EXPECT_TRUE(typeid(vector_field_test) == typeid(vector_field));
}


// Tests on DefaultHostExecutionSpace ------------------------------------------------------------
TEST(MemorySpace, FieldMemOnHostT)
{
    host_t<FieldMemX<double>> field_test_alloc(idx_range_x);
    DFieldMem<IdxRangeX, ddc::KokkosAllocator<double, Kokkos::HostSpace>> field_alloc(idx_range_x);
    EXPECT_TRUE(typeid(field_test_alloc) == typeid(field_alloc));
}


TEST(MemorySpace, FieldOnHostT)
{
    host_t<FieldMemX<double>> field_test_alloc(idx_range_x);
    host_t<FieldX<double>> field_test = get_field(field_test_alloc);
    host_t<DFieldMem<IdxRangeX, ddc::KokkosAllocator<double, Kokkos::HostSpace>>> field_alloc(
            idx_range_x);
    host_t<DField<IdxRangeX>> field = get_field(field_alloc);
    EXPECT_TRUE(typeid(field_test) == typeid(field));
}


TEST(MemorySpace, VectorFieldMemOnHostT)
{
    host_t<VectorFieldMemX<double>> vector_field_test_alloc(idx_range_x);
    VectorFieldMem<double, IdxRangeX, Direction, ddc::KokkosAllocator<double, Kokkos::HostSpace>>
            vector_field_alloc(idx_range_x);
    EXPECT_TRUE(typeid(vector_field_test_alloc) == typeid(vector_field_alloc));
}


TEST(MemorySpace, VectorFieldOnHostT)
{
    host_t<VectorFieldMemX<double>> vector_field_test_alloc(idx_range_x);
    host_t<VectorFieldX<double>> vector_field_test = get_field(vector_field_test_alloc);
    VectorFieldMem<double, IdxRangeX, Direction, ddc::KokkosAllocator<double, Kokkos::HostSpace>>
            vector_field_alloc(idx_range_x);
    VectorField vector_field = get_field(vector_field_alloc);
    EXPECT_TRUE(typeid(vector_field_test) == typeid(vector_field));
}
