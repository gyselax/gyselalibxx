#include <algorithm>
#include <cstdint>
#include <cstring>
#include <iosfwd>
#include <vector>

#include <benchmark/benchmark.h>

#include "block_span.h"
#include "deepcopy.h"
#include "mcoord.h"
#include "mdomain.h"
#include "product_mdomain.h"
#include "product_mesh.h"
#include "rcoord.h"
#include "taggedvector.h"
#include "uniform_mesh.h"

class DimX;
class DimVx;

using MeshX = UniformMesh<DimX>;
using MDomainX = ProductMDomain<MeshX>;
using DBlockX = Block<MDomainX, double>;
using DBlockSpanX = BlockSpan<MDomainX, double>;
using MCoordX = MCoord<MeshX>;
using RCoordX = RCoord<DimX>;

using MeshVx = UniformMesh<DimVx>;
using MDomainVx = ProductMDomain<MeshVx>;
using DBlockVx = Block<MDomainVx, double>;
using DBlockSpanVx = BlockSpan<MDomainVx, double>;
using MCoordVx = MCoord<MeshVx>;
using RCoordVx = RCoord<DimVx>;

using MeshXVx = ProductMesh<MeshX, MeshVx>;
using MDomainXVx = ProductMDomain<MeshX, MeshVx>;
using DBlockXVx = Block<MDomainXVx, double>;
using DBlockSpanXVx = BlockSpan<MDomainXVx, double>;
using MCoordXVx = MCoord<MeshX, MeshVx>;
using RCoordXVx = RCoord<DimX, DimVx>;

using MeshVxX = ProductMesh<MeshVx, MeshX>;
using MDomainVxX = ProductMDomain<MeshVx, MeshX>;
using DBlockVxX = Block<MDomainVxX, double>;
using DBlockSpanVxX = BlockSpan<MDomainVxX, double>;
using MCoordVxX = MCoord<MeshVx, MeshX>;
using RCoordVxX = RCoord<DimVx, DimX>;

static void memcpy_1d(benchmark::State& state)
{
    std::vector<double> src_data(state.range(0), 0.0);
    std::vector<double> dst_data(state.range(0), -1.0);
    for (auto _ : state) {
        std::memcpy(dst_data.data(), src_data.data(), dst_data.size() * sizeof(double));
    }
    state.SetBytesProcessed(int64_t(state.iterations()) * int64_t(state.range(0) * sizeof(double)));
}

static void deepcopy_1d(benchmark::State& state)
{
    std::vector<double> src_data(state.range(0), 0.0);
    std::vector<double> dst_data(state.range(0), -1.0);
    MeshVx mesh_vx(RCoordVx(0.), RCoordVx(2.), state.range(0));
    ProductMesh mesh(mesh_vx);
    MDomainVx const dom(mesh, MCoordVx(state.range(0)));
    DBlockSpanVx src(dom, src_data.data());
    DBlockSpanVx dst(dom, dst_data.data());
    for (auto _ : state) {
        deepcopy(dst, src);
    }
    state.SetBytesProcessed(int64_t(state.iterations()) * int64_t(state.range(0) * sizeof(double)));
}

static void memcpy_2d(benchmark::State& state)
{
    std::vector<double> src_data(state.range(0) * state.range(1), 0.0);
    std::vector<double> dst_data(state.range(0) * state.range(1), -1.0);
    for (auto _ : state) {
        for (std::size_t i = 0; i < state.range(0); ++i) {
            std::
                    memcpy(dst_data.data() + i * state.range(1),
                           src_data.data() + i * state.range(1),
                           state.range(1) * sizeof(double));
        }
    }
    state.SetBytesProcessed(
            int64_t(state.iterations())
            * int64_t(state.range(0) * state.range(1) * sizeof(double)));
}

static void deepcopy_2d(benchmark::State& state)
{
    std::vector<double> src_data(state.range(0) * state.range(1), 0.0);
    std::vector<double> dst_data(state.range(0) * state.range(1), -1.0);
    MeshX mesh_x(RCoordX(0.), RCoordX(2.), state.range(0));
    MeshVx mesh_vx(RCoordVx(0.), RCoordVx(2.), state.range(1));
    ProductMesh mesh(mesh_x, mesh_vx);
    MDomainXVx const dom(mesh, MCoordXVx(state.range(0) - 1, state.range(1) - 1));
    DBlockSpanXVx src(dom, src_data.data());
    DBlockSpanXVx dst(dom, dst_data.data());
    for (auto _ : state) {
        deepcopy(dst, src);
    }
    state.SetBytesProcessed(
            int64_t(state.iterations())
            * int64_t(state.range(0) * state.range(1) * sizeof(double)));
}

static void deepcopy_subblock_2d(benchmark::State& state)
{
    std::vector<double> src_data(state.range(0) * state.range(1), 0.0);
    std::vector<double> dst_data(state.range(0) * state.range(1), -1.0);
    MeshX mesh_x(RCoordX(0.), RCoordX(2.), state.range(0));
    MeshVx mesh_vx(RCoordVx(0.), RCoordVx(2.), state.range(1));
    ProductMesh mesh(mesh_x, mesh_vx);
    MDomainXVx const dom(mesh, MCoordXVx(state.range(0) - 1, state.range(1) - 1));
    DBlockSpanXVx src(dom, src_data.data());
    DBlockSpanXVx dst(dom, dst_data.data());
    for (auto _ : state) {
        for (MCoordX i : get<MeshX>(dom)) {
            auto&& dst_i = dst[i];
            auto&& src_i = src[i];
            deepcopy(dst_i, src_i);
        }
    }
    state.SetBytesProcessed(
            int64_t(state.iterations())
            * int64_t(state.range(0) * state.range(1) * sizeof(double)));
}

// Let say 1MB cache
constexpr std::size_t small_dim1_2D = 400;
constexpr std::size_t small_dim2_2D = small_dim1_2D;

constexpr std::size_t small_dim1_1D = small_dim1_2D * small_dim1_2D;

constexpr std::size_t large_dim1_2D = 2000;
constexpr std::size_t large_dim2_2D = large_dim1_2D;

constexpr std::size_t large_dim1_1D = large_dim1_2D * large_dim1_2D;

// 1D
BENCHMARK(memcpy_1d)->Arg(small_dim1_1D);
BENCHMARK(deepcopy_1d)->Arg(small_dim1_1D);
BENCHMARK(memcpy_1d)->Arg(large_dim1_1D);
BENCHMARK(deepcopy_1d)->Arg(large_dim1_1D);

// 2D
BENCHMARK(memcpy_2d)->Args({small_dim1_2D, small_dim2_2D});
BENCHMARK(deepcopy_subblock_2d)->Args({small_dim1_2D, small_dim2_2D});
BENCHMARK(memcpy_2d)->Args({large_dim1_2D, large_dim2_2D});
BENCHMARK(deepcopy_2d)->Args({large_dim1_2D, large_dim2_2D});
BENCHMARK(deepcopy_subblock_2d)->Args({large_dim1_2D, large_dim2_2D});

BENCHMARK_MAIN();
