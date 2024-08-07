// SPDX-License-Identifier: MIT
#include <typeinfo>

#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

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


struct DDimX
{
};
using DElemX = ddc::DiscreteElement<DDimX>;
using DVectX = ddc::DiscreteVector<DDimX>;
using DDomX = ddc::DiscreteDomain<DDimX>;


template <class Datatype>
using ChunkX = ddc::Chunk<Datatype, DDomX>;
template <class Datatype>
using ChunkSpanX = ddc::ChunkSpan<Datatype, DDomX>;


template <class Datatype>
using VectorFieldX = VectorField<Datatype, DDomX, Direction>;
template <class Datatype>
using VectorFieldSpanX = VectorFieldSpan<Datatype, DDomX, Direction>;


static DElemX constexpr lbound_x(50);
static DVectX constexpr nelems_x(3);
static DDomX constexpr dom_x(lbound_x, nelems_x);

} // end namespace


// Tests on DefaultExecutionSpace ----------------------------------------------------------------
TEST(MemorySpace, ChunkOnDeviceT)
{
    device_t<ChunkX<double>> chunk_test_alloc(dom_x);
    ddc::Chunk<
            double,
            DDomX,
            ddc::KokkosAllocator<double, Kokkos::DefaultExecutionSpace::memory_space>>
            chunk_alloc(dom_x);
    EXPECT_TRUE(typeid(chunk_test_alloc) == typeid(chunk_alloc));
}


TEST(MemorySpace, ChunkSpanOnDeviceT)
{
    device_t<ChunkX<double>> chunk_test_alloc(dom_x);
    device_t<ChunkSpanX<double>> chunk_span_test = chunk_test_alloc.span_view();
    ddc::Chunk<
            double,
            DDomX,
            ddc::KokkosAllocator<double, Kokkos::DefaultExecutionSpace::memory_space>>
            chunk_alloc(dom_x);
    ddc::ChunkSpan chunk_span = chunk_alloc.span_view();
    EXPECT_TRUE(typeid(chunk_span_test) == typeid(chunk_span));
}


TEST(MemorySpace, VectorFieldOnDeviceT)
{
    device_t<VectorFieldX<double>> vector_field_test_alloc(dom_x);
    VectorField<
            double,
            DDomX,
            Direction,
            ddc::KokkosAllocator<double, Kokkos::DefaultExecutionSpace::memory_space>>
            vector_field_alloc(dom_x);
    EXPECT_TRUE(typeid(vector_field_test_alloc) == typeid(vector_field_alloc));
}


TEST(MemorySpace, VectorFieldSpanOnDeviceT)
{
    device_t<VectorFieldX<double>> vector_field_test_alloc(dom_x);
    device_t<VectorFieldSpanX<double>> vector_field_span_test = vector_field_test_alloc.span_view();
    VectorField<
            double,
            DDomX,
            Direction,
            ddc::KokkosAllocator<double, Kokkos::DefaultExecutionSpace::memory_space>>
            vector_field_alloc(dom_x);
    VectorFieldSpan vector_field_span = vector_field_alloc.span_view();
    EXPECT_TRUE(typeid(vector_field_span_test) == typeid(vector_field_span));
}


// Tests on DefaultHostExecutionSpace ------------------------------------------------------------
TEST(MemorySpace, ChunkOnHostT)
{
    host_t<ChunkX<double>> chunk_test_alloc(dom_x);
    ddc::Chunk<
            double,
            DDomX,
            ddc::KokkosAllocator<double, Kokkos::DefaultHostExecutionSpace::memory_space>>
            chunk_alloc(dom_x);
    EXPECT_TRUE(typeid(chunk_test_alloc) == typeid(chunk_alloc));
}


TEST(MemorySpace, ChunkSpanOnHostT)
{
    host_t<ChunkX<double>> chunk_test_alloc(dom_x);
    host_t<ChunkSpanX<double>> chunk_span_test = chunk_test_alloc.span_view();
    ddc::Chunk<
            double,
            DDomX,
            ddc::KokkosAllocator<double, Kokkos::DefaultHostExecutionSpace::memory_space>>
            chunk_alloc(dom_x);
    ddc::ChunkSpan chunk_span = chunk_alloc.span_view();
    EXPECT_TRUE(typeid(chunk_span_test) == typeid(chunk_span));
}


TEST(MemorySpace, VectorFieldOnHostT)
{
    host_t<VectorFieldX<double>> vector_field_test_alloc(dom_x);
    VectorField<
            double,
            DDomX,
            Direction,
            ddc::KokkosAllocator<double, Kokkos::DefaultHostExecutionSpace::memory_space>>
            vector_field_alloc(dom_x);
    EXPECT_TRUE(typeid(vector_field_test_alloc) == typeid(vector_field_alloc));
}


TEST(MemorySpace, VectorFieldSpanOnHostT)
{
    host_t<VectorFieldX<double>> vector_field_test_alloc(dom_x);
    host_t<VectorFieldSpanX<double>> vector_field_span_test = vector_field_test_alloc.span_view();
    VectorField<
            double,
            DDomX,
            Direction,
            ddc::KokkosAllocator<double, Kokkos::DefaultHostExecutionSpace::memory_space>>
            vector_field_alloc(dom_x);
    VectorFieldSpan vector_field_span = vector_field_alloc.span_view();
    EXPECT_TRUE(typeid(vector_field_span_test) == typeid(vector_field_span));
}