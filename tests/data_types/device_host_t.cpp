// SPDX-License-Identifier: MIT
#include <typeinfo>

#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
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
using VectorFieldX = VectorField<Datatype, IdxRangeX, Direction>;
template <class Datatype>
using VectorFieldFieldX = VectorFieldSpan<Datatype, IdxRangeX, Direction>;


static IdxX constexpr lbound_x(50);
static IdxSteptX constexpr nelems_x(3);
static IdxRangeX constexpr dom_x(lbound_x, nelems_x);

} // end namespace


// Tests on DefaultExecutionSpace ----------------------------------------------------------------
TEST(MemorySpace, ChunkOnDeviceT)
{
    FieldMemX<double> chunk_test_alloc(dom_x);
    FieldMem<
            double,
            IdxRangeX,
            ddc::KokkosAllocator<double, Kokkos::DefaultExecutionSpace::memory_space>>
            chunk_alloc(dom_x);
    EXPECT_TRUE(typeid(chunk_test_alloc) == typeid(chunk_alloc));
}


TEST(MemorySpace, ChunkSpanOnDeviceT)
{
    FieldMemX<double> chunk_test_alloc(dom_x);
    FieldX<double> chunk_span_test = get_field(chunk_test_alloc);
    FieldMem<
            double,
            IdxRangeX,
            ddc::KokkosAllocator<double, Kokkos::DefaultExecutionSpace::memory_space>>
            chunk_alloc(dom_x);
    ddc::ChunkSpan chunk_span = get_field(chunk_alloc);
    EXPECT_TRUE(typeid(chunk_span_test) == typeid(chunk_span));
}


TEST(MemorySpace, VectorFieldOnDeviceT)
{
    device_t<VectorFieldX<double>> vector_field_test_alloc(dom_x);
    VectorField<
            double,
            IdxRangeX,
            Direction,
            ddc::KokkosAllocator<double, Kokkos::DefaultExecutionSpace::memory_space>>
            vector_field_alloc(dom_x);
    EXPECT_TRUE(typeid(vector_field_test_alloc) == typeid(vector_field_alloc));
}


TEST(MemorySpace, VectorFieldSpanOnDeviceT)
{
    device_t<VectorFieldX<double>> vector_field_test_alloc(dom_x);
    device_t<VectorFieldFieldX<double>> vector_field_span_test = get_field(vector_field_test_alloc);
    VectorField<
            double,
            IdxRangeX,
            Direction,
            ddc::KokkosAllocator<double, Kokkos::DefaultExecutionSpace::memory_space>>
            vector_field_alloc(dom_x);
    VectorFieldSpan vector_field_span = get_field(vector_field_alloc);
    EXPECT_TRUE(typeid(vector_field_span_test) == typeid(vector_field_span));
}


// Tests on DefaultHostExecutionSpace ------------------------------------------------------------
TEST(MemorySpace, ChunkOnHostT)
{
    host_t<FieldMemX<double>> chunk_test_alloc(dom_x);
    FieldMem<
            double,
            IdxRangeX,
            ddc::KokkosAllocator<double, Kokkos::DefaultHostExecutionSpace::memory_space>>
            chunk_alloc(dom_x);
    EXPECT_TRUE(typeid(chunk_test_alloc) == typeid(chunk_alloc));
}


TEST(MemorySpace, ChunkSpanOnHostT)
{
    host_t<FieldMemX<double>> chunk_test_alloc(dom_x);
    host_t<FieldX<double>> chunk_span_test = get_field(chunk_test_alloc);
    FieldMem<
            double,
            IdxRangeX,
            ddc::KokkosAllocator<double, Kokkos::DefaultHostExecutionSpace::memory_space>>
            chunk_alloc(dom_x);
    ddc::ChunkSpan chunk_span = get_field(chunk_alloc);
    EXPECT_TRUE(typeid(chunk_span_test) == typeid(chunk_span));
}


TEST(MemorySpace, VectorFieldOnHostT)
{
    host_t<VectorFieldX<double>> vector_field_test_alloc(dom_x);
    VectorField<
            double,
            IdxRangeX,
            Direction,
            ddc::KokkosAllocator<double, Kokkos::DefaultHostExecutionSpace::memory_space>>
            vector_field_alloc(dom_x);
    EXPECT_TRUE(typeid(vector_field_test_alloc) == typeid(vector_field_alloc));
}


TEST(MemorySpace, VectorFieldSpanOnHostT)
{
    host_t<VectorFieldX<double>> vector_field_test_alloc(dom_x);
    host_t<VectorFieldFieldX<double>> vector_field_span_test = get_field(vector_field_test_alloc);
    VectorField<
            double,
            IdxRangeX,
            Direction,
            ddc::KokkosAllocator<double, Kokkos::DefaultHostExecutionSpace::memory_space>>
            vector_field_alloc(dom_x);
    VectorFieldSpan vector_field_span = get_field(vector_field_alloc);
    EXPECT_TRUE(typeid(vector_field_span_test) == typeid(vector_field_span));
}
