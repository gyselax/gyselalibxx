// SPDX-License-Identifier: MIT
#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include <mpilayout.hpp>
#include <mpitransposealltoall.hpp>

#include "ddc_alias_inline_functions.hpp"

namespace {

struct W
{
};
struct X
{
};
struct Y
{
};
struct Z
{
};

struct GridW : UniformGridBase<W>
{
};
struct GridX : UniformGridBase<X>
{
};
struct GridY : UniformGridBase<Y>
{
};
struct GridZ : UniformGridBase<Z>
{
};

using IdxW = Idx<GridW>;
using IdxX = Idx<GridX>;
using IdxY = Idx<GridY>;
using IdxZ = Idx<GridZ>;
using IdxXY = Idx<GridX, GridY>;
using IdxYX = Idx<GridY, GridX>;
using IdxXYZ = Idx<GridX, GridY, GridZ>;
using IdxYZX = Idx<GridY, GridZ, GridX>;
using IdxWXYZ = Idx<GridW, GridX, GridY, GridZ>;
using IdxYZWX = Idx<GridY, GridZ, GridW, GridX>;

using IdxStepW = IdxStep<GridW>;
using IdxStepX = IdxStep<GridX>;
using IdxStepY = IdxStep<GridY>;
using IdxStepZ = IdxStep<GridZ>;
using IdxStepXY = IdxStep<GridX, GridY>;
using IdxStepXYZ = IdxStep<GridX, GridY, GridZ>;
using IdxStepWXYZ = IdxStep<GridW, GridX, GridY, GridZ>;

using IdxRangeW = IdxRange<GridW>;
using IdxRangeX = IdxRange<GridX>;
using IdxRangeY = IdxRange<GridY>;
using IdxRangeZ = IdxRange<GridZ>;
using IdxRangeXY = IdxRange<GridX, GridY>;
using IdxRangeYX = IdxRange<GridY, GridX>;
using IdxRangeXYZ = IdxRange<GridX, GridY, GridZ>;
using IdxRangeYZX = IdxRange<GridY, GridZ, GridX>;
using IdxRangeWXYZ = IdxRange<GridW, GridX, GridY, GridZ>;
using IdxRangeYZWX = IdxRange<GridY, GridZ, GridW, GridX>;

using CoordW = Coord<W>;
using CoordX = Coord<X>;
using CoordY = Coord<Y>;
using CoordZ = Coord<Z>;

using IFieldMemXY = host_t<FieldMem<std::size_t, IdxRangeXY>>;
using IFieldMemYX = host_t<FieldMem<std::size_t, IdxRangeYX>>;
using IFieldMemXYZ = host_t<FieldMem<std::size_t, IdxRangeXYZ>>;
using IFieldMemYZX = host_t<FieldMem<std::size_t, IdxRangeYZX>>;
using IFieldMemWXYZ = host_t<FieldMem<std::size_t, IdxRangeWXYZ>>;
using IFieldMemYZWX = host_t<FieldMem<std::size_t, IdxRangeYZWX>>;

using IFieldYX = host_t<Field<std::size_t, IdxRangeYX>>;

using XDistribLayout = MPILayout<IdxRangeXY, GridX>;
using YDistribLayout = MPILayout<IdxRangeYX, GridY>;

using YDistribLayout3D = MPILayout<IdxRangeXYZ, GridY>;
using ZDistribLayout3D = MPILayout<IdxRangeYZX, GridZ>;

using XYDistribLayout4D = MPILayout<IdxRangeWXYZ, GridX, GridY>;
using WZDistribLayout4D = MPILayout<IdxRangeYZWX, GridW, GridZ>;

template <class HeadDim, class... Dims>
KOKKOS_FUNCTION std::size_t get_unique_id(
        Idx<HeadDim, Dims...> index,
        IdxRange<HeadDim, Dims...> global_size)
{
    if constexpr (sizeof...(Dims) == 0) {
        return (index - Idx<HeadDim>(0)).value();
    } else {
        Idx<HeadDim> head_index(index);
        IdxRange<Dims...> lower_sizes(global_size);
        Idx<Dims...> lower_indices(index);
        std::size_t head_index_val = (head_index - Idx<HeadDim>(0)).value();
        return head_index_val * lower_sizes.size() + get_unique_id(lower_indices, lower_sizes);
    }
}

void test_AllToAll2D_GPU()
{
    IdxStepX x_size(10);
    IdxStepY y_size(12);

    IdxXY idx_range_start(0, 0);
    IdxStepXY dom_size(x_size, y_size);
    IdxRangeXY full_idx_range(idx_range_start, dom_size);

    MPITransposeAllToAll<XDistribLayout, YDistribLayout> transpose(full_idx_range, MPI_COMM_WORLD);

    host_t<IFieldMemXY> recv_buffer_host(transpose.get_local_idx_range<XDistribLayout>());
    device_t<IFieldMemXY> recv_buffer(transpose.get_local_idx_range<XDistribLayout>());
    device_t<IFieldMemYX> send_buffer_alloc(transpose.get_local_idx_range<YDistribLayout>());
    device_t<IFieldYX> send_buffer = get_field(send_buffer_alloc);

    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            get_idx_range(send_buffer),
            KOKKOS_LAMBDA(IdxYX ixy) {
                send_buffer(ixy) = get_unique_id(IdxXY(ixy), full_idx_range);
            });

    transpose(
            Kokkos::DefaultExecutionSpace(),
            get_field(recv_buffer),
            get_const_field(send_buffer));

    ddc::parallel_deepcopy(get_field(recv_buffer_host), get_const_field(recv_buffer));

    bool success = true;
    ddc::for_each(get_idx_range(recv_buffer_host), [&](IdxXY ixy) {
        success = success and (recv_buffer_host(ixy) == get_unique_id(ixy, full_idx_range));
    });
    EXPECT_TRUE(success);
}

} // namespace

TEST(MPIParallelisation, AllToAll2D_CPU)
{
    IdxStepX x_size(10);
    IdxStepY y_size(12);

    IdxXY idx_range_start(0, 0);
    IdxStepXY dom_size(x_size, y_size);
    IdxRangeXY full_idx_range(idx_range_start, dom_size);

    MPITransposeAllToAll<XDistribLayout, YDistribLayout> transpose(full_idx_range, MPI_COMM_WORLD);

    IFieldMemXY recv_buffer(transpose.get_local_idx_range<XDistribLayout>());
    IFieldMemYX send_buffer(transpose.get_local_idx_range<YDistribLayout>());

    ddc::for_each(get_idx_range(send_buffer), [&](IdxYX ixy) {
        send_buffer(ixy) = get_unique_id(IdxXY(ixy), full_idx_range);
    });

    transpose(
            Kokkos::DefaultHostExecutionSpace(),
            get_field(recv_buffer),
            get_const_field(send_buffer));

    bool success = true;
    ddc::for_each(get_idx_range(recv_buffer), [&](IdxXY ixy) {
        success = success and (recv_buffer(ixy) == get_unique_id(ixy, full_idx_range));
    });
    EXPECT_TRUE(success);
}

TEST(MPIParallelisation, AllToAll2D_GPU)
{
    test_AllToAll2D_GPU();
}

TEST(MPIParallelisation, AllToAll3D_CPU)
{
    IdxStepX x_size(10);
    IdxStepY y_size(12);
    IdxStepZ z_size(4);

    IdxXYZ idx_range_start(0, 0, 0);
    IdxStepXYZ dom_size(x_size, y_size, z_size);
    IdxRangeXYZ full_idx_range(idx_range_start, dom_size);

    MPITransposeAllToAll<YDistribLayout3D, ZDistribLayout3D>
            transpose(full_idx_range, MPI_COMM_WORLD);

    IFieldMemXYZ send_buffer(transpose.get_local_idx_range<YDistribLayout3D>());
    IFieldMemYZX recv_buffer(transpose.get_local_idx_range<ZDistribLayout3D>());

    ddc::for_each(get_idx_range(send_buffer), [&](IdxXYZ ixyz) {
        send_buffer(ixyz) = get_unique_id(ixyz, full_idx_range);
    });

    MPI_Barrier(MPI_COMM_WORLD);

    transpose(
            Kokkos::DefaultHostExecutionSpace(),
            get_field(recv_buffer),
            get_const_field(send_buffer));

    bool success = true;
    ddc::for_each(get_idx_range(recv_buffer), [&](IdxYZX ixyz) {
        double expected = get_unique_id(IdxXYZ(ixyz), full_idx_range);
        success = success and (recv_buffer(ixyz) == expected);
    });
    EXPECT_TRUE(success);
}

TEST(MPIParallelisation, AllToAll4D_CPU)
{
    IdxStepW w_size(10);
    IdxStepX x_size(10);
    IdxStepY y_size(10);
    IdxStepZ z_size(10);

    IdxWXYZ idx_range_start(0, 0, 0, 0);
    IdxStepWXYZ dom_size(w_size, x_size, y_size, z_size);
    IdxRangeWXYZ full_idx_range(idx_range_start, dom_size);

    MPITransposeAllToAll<XYDistribLayout4D, WZDistribLayout4D>
            transpose(full_idx_range, MPI_COMM_WORLD);

    IFieldMemWXYZ send_buffer(transpose.get_local_idx_range<XYDistribLayout4D>());
    IFieldMemYZWX recv_buffer(transpose.get_local_idx_range<WZDistribLayout4D>());

    ddc::for_each(get_idx_range(send_buffer), [&](IdxWXYZ iwxyz) {
        send_buffer(iwxyz) = get_unique_id(iwxyz, full_idx_range);
    });

    MPI_Barrier(MPI_COMM_WORLD);

    transpose(
            Kokkos::DefaultHostExecutionSpace(),
            get_field(recv_buffer),
            get_const_field(send_buffer));

    bool success = true;
    ddc::for_each(get_idx_range(recv_buffer), [&](IdxYZWX iwxyz) {
        std::size_t expected = get_unique_id(IdxWXYZ(iwxyz), full_idx_range);
        success = success and (recv_buffer(iwxyz) == expected);
    });
    EXPECT_TRUE(success);
}
