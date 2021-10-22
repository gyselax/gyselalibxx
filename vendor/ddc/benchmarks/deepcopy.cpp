// SPDX-License-Identifier: MIT
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <iosfwd>
#include <vector>

#include <ddc/ChunkSpan>
#include <ddc/Coordinate>
#include <ddc/DiscreteCoordinate>
#include <ddc/DiscreteDomain>
#include <ddc/UniformDiscretization>

#include <benchmark/benchmark.h>

namespace {

class DimX;
class DimY;

using DDimX = UniformDiscretization<DimX>;
using IDomainX = DiscreteDomain<DDimX>;
using DChunkX = Chunk<double, IDomainX>;
using DChunkSpanX = ChunkSpan<double, IDomainX>;
using IndexX = DiscreteCoordinate<DDimX>;
using CoordX = Coordinate<DimX>;

using DDimY = UniformDiscretization<DimY>;
using MDomainY = DiscreteDomain<DDimY>;
using DChunkY = Chunk<double, MDomainY>;
using DChunkSpanY = ChunkSpan<double, MDomainY>;
using MCoordY = DiscreteCoordinate<DDimY>;
using MLengthY = DiscreteVector<DDimY>;
using RCoordY = Coordinate<DimY>;

using MDomainSpXY = DiscreteDomain<DDimX, DDimY>;
using DChunkSpXY = Chunk<double, MDomainSpXY>;
using DChunkSpanXY = ChunkSpan<double, MDomainSpXY>;
using MCoordXY = DiscreteCoordinate<DDimX, DDimY>;
using MLengthXY = DiscreteVector<DDimX, DDimY>;
using RCoordXY = Coordinate<DimX, DimY>;

using MDomainYX = DiscreteDomain<DDimY, DDimX>;
using DChunkYX = Chunk<double, MDomainYX>;
using DChunkSpanYX = ChunkSpan<double, MDomainYX>;
using MCoordYX = DiscreteCoordinate<DDimY, DDimX>;
using RCoordYX = Coordinate<DimY, DimX>;

} // namespace

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
    DDimY ddim_y(RCoordY(0.), RCoordY(2.), state.range(0));
    MDomainY const dom(ddim_y, MLengthY(state.range(0)));
    DChunkSpanY src(src_data.data(), dom);
    DChunkSpanY dst(dst_data.data(), dom);
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
    DDimX ddim_x(CoordX(0.), CoordX(2.), state.range(0));
    DDimY ddim_y(RCoordY(0.), RCoordY(2.), state.range(1));
    MDomainSpXY const dom(ddim_x, ddim_y, MLengthXY(state.range(0) - 1, state.range(1) - 1));
    DChunkSpanXY src(src_data.data(), dom);
    DChunkSpanXY dst(dst_data.data(), dom);
    for (auto _ : state) {
        deepcopy(dst, src);
    }
    state.SetBytesProcessed(
            int64_t(state.iterations())
            * int64_t(state.range(0) * state.range(1) * sizeof(double)));
}

static void deepcopy_subchunck_2d(benchmark::State& state)
{
    std::vector<double> src_data(state.range(0) * state.range(1), 0.0);
    std::vector<double> dst_data(state.range(0) * state.range(1), -1.0);
    DDimX ddim_x(CoordX(0.), CoordX(2.), state.range(0));
    DDimY ddim_y(RCoordY(0.), RCoordY(2.), state.range(1));
    MDomainSpXY const dom(ddim_x, ddim_y, MLengthXY(state.range(0) - 1, state.range(1) - 1));
    DChunkSpanXY src(src_data.data(), dom);
    DChunkSpanXY dst(dst_data.data(), dom);
    for (auto _ : state) {
        for (IndexX i : select<DDimX>(dom)) {
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
std::size_t constexpr small_dim1_2D = 400;
std::size_t constexpr small_dim2_2D = small_dim1_2D;

std::size_t constexpr small_dim1_1D = small_dim1_2D * small_dim1_2D;

std::size_t constexpr large_dim1_2D = 2000;
std::size_t constexpr large_dim2_2D = large_dim1_2D;

std::size_t constexpr large_dim1_1D = large_dim1_2D * large_dim1_2D;

// 1D
BENCHMARK(memcpy_1d)->Arg(small_dim1_1D);
BENCHMARK(deepcopy_1d)->Arg(small_dim1_1D);
BENCHMARK(memcpy_1d)->Arg(large_dim1_1D);
BENCHMARK(deepcopy_1d)->Arg(large_dim1_1D);

// 2D
BENCHMARK(memcpy_2d)->Args({small_dim1_2D, small_dim2_2D});
BENCHMARK(deepcopy_subchunck_2d)->Args({small_dim1_2D, small_dim2_2D});
BENCHMARK(memcpy_2d)->Args({large_dim1_2D, large_dim2_2D});
BENCHMARK(deepcopy_2d)->Args({large_dim1_2D, large_dim2_2D});
BENCHMARK(deepcopy_subchunck_2d)->Args({large_dim1_2D, large_dim2_2D});

BENCHMARK_MAIN();
