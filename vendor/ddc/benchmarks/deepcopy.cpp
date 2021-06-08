#include <cstring>
#include <vector>

#include <benchmark/benchmark.h>

#include "blockview.h"

class DimX;
class DimVx;

using MeshX = UniformMesh<DimX>;
using MDomainX = UniformMDomain<DimX>;
using DBlockX = Block<MDomainX, double>;
using DBlockSpanX = BlockView<MDomainX, double>;
using MCoordX = MCoord<DimX>;
using RCoordX = RCoord<DimX>;

using MeshVx = UniformMesh<DimVx>;
using MDomainVx = UniformMDomain<DimVx>;
using DBlockVx = Block<MDomainVx, double>;
using DBlockSpanVx = BlockView<MDomainVx, double>;
using MCoordVx = MCoord<DimVx>;
using RCoordVx = RCoord<DimVx>;

using MeshXVx = UniformMesh<DimX, DimVx>;
using MDomainXVx = UniformMDomain<DimX, DimVx>;
using DBlockXVx = Block<MDomainXVx, double>;
using DBlockSpanXVx = BlockView<MDomainXVx, double>;
using MCoordXVx = MCoord<DimX, DimVx>;
using RCoordXVx = RCoord<DimX, DimVx>;

using MeshVxX = UniformMesh<DimVx, DimX>;
using MDomainVxX = UniformMDomain<DimVx, DimX>;
using DBlockVxX = Block<MDomainVxX, double>;
using DBlockSpanVxX = BlockView<MDomainVxX, double>;
using MCoordVxX = MCoord<DimVx, DimX>;
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
    MDomainVx const dom(RCoordVx(0.), RCoordVx(2.), MCoordVx(0), MCoordVx(state.range(0)));
    DBlockSpanVx src(dom.mesh(), DSpan1D(src_data.data(), state.range(0)));
    DBlockSpanVx dst(dom.mesh(), DSpan1D(dst_data.data(), state.range(0)));
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
    MDomainXVx const
            dom(RCoordXVx(0., 0.),
                RCoordXVx(2., 2.),
                MCoordXVx(0, 0),
                MCoordXVx(state.range(0), state.range(1)));
    DBlockSpanXVx src(dom.mesh(), DSpan2D(src_data.data(), state.range(0), state.range(1)));
    DBlockSpanXVx dst(dom.mesh(), DSpan2D(dst_data.data(), state.range(0), state.range(1)));
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
    MDomainXVx const
            dom(RCoordXVx(0., 0.),
                RCoordXVx(2., 2.),
                MCoordXVx(0, 0),
                MCoordXVx(state.range(0), state.range(1)));
    DBlockSpanXVx src(dom.mesh(), DSpan2D(src_data.data(), state.range(0), state.range(1)));
    DBlockSpanXVx dst(dom.mesh(), DSpan2D(dst_data.data(), state.range(0), state.range(1)));
    for (auto _ : state) {
        for (std::size_t i = 0; i < src.extent(0); ++i) {
            auto&& dst_i = dst.subblockview(i, std::experimental::all);
            auto&& src_i = src.subblockview(i, std::experimental::all);
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
