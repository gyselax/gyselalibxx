// SPDX-License-Identifier: MIT
#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include <mpilayout.hpp>
#include <mpitransposealltoall.hpp>

namespace {

struct RDimW
{
};
struct RDimX
{
};
struct RDimY
{
};
struct RDimZ
{
};

struct IDimW : ddc::UniformPointSampling<RDimW>
{
};
struct IDimX : ddc::UniformPointSampling<RDimX>
{
};
struct IDimY : ddc::UniformPointSampling<RDimY>
{
};
struct IDimZ : ddc::UniformPointSampling<RDimZ>
{
};

using IndexW = ddc::DiscreteElement<IDimW>;
using IndexX = ddc::DiscreteElement<IDimX>;
using IndexY = ddc::DiscreteElement<IDimY>;
using IndexZ = ddc::DiscreteElement<IDimZ>;
using IndexXY = ddc::DiscreteElement<IDimX, IDimY>;
using IndexYX = ddc::DiscreteElement<IDimY, IDimX>;
using IndexXYZ = ddc::DiscreteElement<IDimX, IDimY, IDimZ>;
using IndexYZX = ddc::DiscreteElement<IDimY, IDimZ, IDimX>;
using IndexWXYZ = ddc::DiscreteElement<IDimW, IDimX, IDimY, IDimZ>;
using IndexYZWX = ddc::DiscreteElement<IDimY, IDimZ, IDimW, IDimX>;

using IVectW = ddc::DiscreteVector<IDimW>;
using IVectX = ddc::DiscreteVector<IDimX>;
using IVectY = ddc::DiscreteVector<IDimY>;
using IVectZ = ddc::DiscreteVector<IDimZ>;
using IVectXY = ddc::DiscreteVector<IDimX, IDimY>;
using IVectXYZ = ddc::DiscreteVector<IDimX, IDimY, IDimZ>;
using IVectWXYZ = ddc::DiscreteVector<IDimW, IDimX, IDimY, IDimZ>;

using IDomainW = ddc::DiscreteDomain<IDimW>;
using IDomainX = ddc::DiscreteDomain<IDimX>;
using IDomainY = ddc::DiscreteDomain<IDimY>;
using IDomainZ = ddc::DiscreteDomain<IDimZ>;
using IDomainXY = ddc::DiscreteDomain<IDimX, IDimY>;
using IDomainYX = ddc::DiscreteDomain<IDimY, IDimX>;
using IDomainXYZ = ddc::DiscreteDomain<IDimX, IDimY, IDimZ>;
using IDomainYZX = ddc::DiscreteDomain<IDimY, IDimZ, IDimX>;
using IDomainWXYZ = ddc::DiscreteDomain<IDimW, IDimX, IDimY, IDimZ>;
using IDomainYZWX = ddc::DiscreteDomain<IDimY, IDimZ, IDimW, IDimX>;

using CoordW = ddc::Coordinate<RDimW>;
using CoordX = ddc::Coordinate<RDimX>;
using CoordY = ddc::Coordinate<RDimY>;
using CoordZ = ddc::Coordinate<RDimZ>;

using IFieldXY = ddc::Chunk<std::size_t, IDomainXY>;
using IFieldYX = ddc::Chunk<std::size_t, IDomainYX>;
using IFieldXYZ = ddc::Chunk<std::size_t, IDomainXYZ>;
using IFieldYZX = ddc::Chunk<std::size_t, IDomainYZX>;
using IFieldWXYZ = ddc::Chunk<std::size_t, IDomainWXYZ>;
using IFieldYZWX = ddc::Chunk<std::size_t, IDomainYZWX>;

using ISpanYX = ddc::ChunkSpan<std::size_t, IDomainYX>;

using XDistribLayout = MPILayout<IDomainXY, IDimX>;
using YDistribLayout = MPILayout<IDomainYX, IDimY>;

using YDistribLayout3D = MPILayout<IDomainXYZ, IDimY>;
using ZDistribLayout3D = MPILayout<IDomainYZX, IDimZ>;

using XYDistribLayout4D = MPILayout<IDomainWXYZ, IDimX, IDimY>;
using WZDistribLayout4D = MPILayout<IDomainYZWX, IDimW, IDimZ>;

template <class HeadDim, class... Dims>
KOKKOS_FUNCTION std::size_t get_unique_id(
        ddc::DiscreteElement<HeadDim, Dims...> index,
        ddc::DiscreteDomain<HeadDim, Dims...> global_size)
{
    if constexpr (sizeof...(Dims) == 0) {
        return (index - ddc::DiscreteElement<HeadDim>(0)).value();
    } else {
        ddc::DiscreteElement<HeadDim> head_index(index);
        ddc::DiscreteDomain<Dims...> lower_sizes(global_size);
        ddc::DiscreteElement<Dims...> lower_indices(index);
        std::size_t head_index_val = (head_index - ddc::DiscreteElement<HeadDim>(0)).value();
        return head_index_val * lower_sizes.size() + get_unique_id(lower_indices, lower_sizes);
    }
}

void test_AllToAll2D_GPU()
{
    IVectX x_size(10);
    IVectY y_size(12);

    IndexXY domain_start(0, 0);
    IVectXY dom_size(x_size, y_size);
    IDomainXY full_domain(domain_start, dom_size);

    MPITransposeAllToAll<XDistribLayout, YDistribLayout> transpose(full_domain, MPI_COMM_WORLD);

    host_t<IFieldXY> recv_buffer_host(transpose.get_local_domain<XDistribLayout>());
    device_t<IFieldXY> recv_buffer(transpose.get_local_domain<XDistribLayout>());
    device_t<IFieldYX> send_buffer_alloc(transpose.get_local_domain<YDistribLayout>());
    device_t<ISpanYX> send_buffer = send_buffer_alloc.span_view();

    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            send_buffer.domain(),
            KOKKOS_LAMBDA(IndexYX ixy) {
                send_buffer(ixy) = get_unique_id(IndexXY(ixy), full_domain);
            });

    transpose(Kokkos::DefaultExecutionSpace(), recv_buffer.span_view(), send_buffer.span_cview());

    ddc::parallel_deepcopy(recv_buffer_host.span_view(), recv_buffer.span_cview());

    bool success = true;
    ddc::for_each(recv_buffer_host.domain(), [&](IndexXY ixy) {
        success = success and (recv_buffer_host(ixy) == get_unique_id(ixy, full_domain));
    });
    EXPECT_TRUE(success);
}

} // namespace

TEST(MPIParallelisation, AllToAll2D_CPU)
{
    IVectX x_size(10);
    IVectY y_size(12);

    IndexXY domain_start(0, 0);
    IVectXY dom_size(x_size, y_size);
    IDomainXY full_domain(domain_start, dom_size);

    MPITransposeAllToAll<XDistribLayout, YDistribLayout> transpose(full_domain, MPI_COMM_WORLD);

    IFieldXY recv_buffer(transpose.get_local_domain<XDistribLayout>());
    IFieldYX send_buffer(transpose.get_local_domain<YDistribLayout>());

    ddc::for_each(send_buffer.domain(), [&](IndexYX ixy) {
        send_buffer(ixy) = get_unique_id(IndexXY(ixy), full_domain);
    });

    transpose(
            Kokkos::DefaultHostExecutionSpace(),
            recv_buffer.span_view(),
            send_buffer.span_cview());

    bool success = true;
    ddc::for_each(recv_buffer.domain(), [&](IndexXY ixy) {
        success = success and (recv_buffer(ixy) == get_unique_id(ixy, full_domain));
    });
    EXPECT_TRUE(success);
}

TEST(MPIParallelisation, AllToAll2D_GPU)
{
    test_AllToAll2D_GPU();
}

TEST(MPIParallelisation, AllToAll3D_CPU)
{
    IVectX x_size(10);
    IVectY y_size(12);
    IVectZ z_size(4);

    IndexXYZ domain_start(0, 0, 0);
    IVectXYZ dom_size(x_size, y_size, z_size);
    IDomainXYZ full_domain(domain_start, dom_size);

    MPITransposeAllToAll<YDistribLayout3D, ZDistribLayout3D> transpose(full_domain, MPI_COMM_WORLD);

    IFieldXYZ send_buffer(transpose.get_local_domain<YDistribLayout3D>());
    IFieldYZX recv_buffer(transpose.get_local_domain<ZDistribLayout3D>());

    ddc::for_each(send_buffer.domain(), [&](IndexXYZ ixyz) {
        send_buffer(ixyz) = get_unique_id(ixyz, full_domain);
    });

    MPI_Barrier(MPI_COMM_WORLD);

    transpose(
            Kokkos::DefaultHostExecutionSpace(),
            recv_buffer.span_view(),
            send_buffer.span_cview());

    bool success = true;
    ddc::for_each(recv_buffer.domain(), [&](IndexYZX ixyz) {
        double expected = get_unique_id(IndexXYZ(ixyz), full_domain);
        success = success and (recv_buffer(ixyz) == expected);
    });
    EXPECT_TRUE(success);
}

TEST(MPIParallelisation, AllToAll4D_CPU)
{
    IVectW w_size(10);
    IVectX x_size(10);
    IVectY y_size(10);
    IVectZ z_size(10);

    IndexWXYZ domain_start(0, 0, 0, 0);
    IVectWXYZ dom_size(w_size, x_size, y_size, z_size);
    IDomainWXYZ full_domain(domain_start, dom_size);

    MPITransposeAllToAll<XYDistribLayout4D, WZDistribLayout4D>
            transpose(full_domain, MPI_COMM_WORLD);

    IFieldWXYZ send_buffer(transpose.get_local_domain<XYDistribLayout4D>());
    IFieldYZWX recv_buffer(transpose.get_local_domain<WZDistribLayout4D>());

    ddc::for_each(send_buffer.domain(), [&](IndexWXYZ iwxyz) {
        send_buffer(iwxyz) = get_unique_id(iwxyz, full_domain);
    });

    MPI_Barrier(MPI_COMM_WORLD);

    transpose(
            Kokkos::DefaultHostExecutionSpace(),
            recv_buffer.span_view(),
            send_buffer.span_cview());

    bool success = true;
    ddc::for_each(recv_buffer.domain(), [&](IndexYZWX iwxyz) {
        std::size_t expected = get_unique_id(IndexWXYZ(iwxyz), full_domain);
        success = success and (recv_buffer(iwxyz) == expected);
    });
    EXPECT_TRUE(success);
}
