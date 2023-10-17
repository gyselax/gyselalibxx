// SPDX-License-Identifier: MIT
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <iosfwd>
#include <vector>

#include <ddc/ddc.hpp>

#include <benchmark/benchmark.h>

namespace {

struct DDimX;
using DElemX = ddc::DiscreteElement<DDimX>;
using DVectX = ddc::DiscreteVector<DDimX>;
using DDomX = ddc::DiscreteDomain<DDimX>;

template <class Datatype>
using ChunkSpanX = ddc::ChunkSpan<Datatype, DDomX>;


struct DDimY;


using DElemXY = ddc::DiscreteElement<DDimX, DDimY>;
using DVectXY = ddc::DiscreteVector<DDimX, DDimY>;
using DDomXY = ddc::DiscreteDomain<DDimX, DDimY>;

template <class Datatype>
using ChunkSpanXY = ddc::ChunkSpan<Datatype, DDomXY>;


// Let say 1MB cache
static std::size_t constexpr small_dim1_2D = 400;
static std::size_t constexpr small_dim2_2D = small_dim1_2D;

static std::size_t constexpr small_dim1_1D = small_dim1_2D * small_dim1_2D;

static std::size_t constexpr large_dim1_2D = 2000;
static std::size_t constexpr large_dim2_2D = large_dim1_2D;

static std::size_t constexpr large_dim1_1D = large_dim1_2D * large_dim1_2D;

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
    DDomX const dom(DElemX(0), DVectX(state.range(0)));
    ChunkSpanX<double> src(src_data.data(), dom);
    ChunkSpanX<double> dst(dst_data.data(), dom);
    for (auto _ : state) {
        ddc::deepcopy(dst, src);
    }
    state.SetBytesProcessed(int64_t(state.iterations()) * int64_t(state.range(0) * sizeof(double)));
}

static void memcpy_2d(benchmark::State& state)
{
    std::vector<double> src_data(state.range(0) * state.range(1), 0.0);
    std::vector<double> dst_data(state.range(0) * state.range(1), -1.0);
    for (auto _ : state) {
        for (int64_t i = 0; i < state.range(0); ++i) {
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
    DDomXY const dom(DElemXY(0, 0), DVectXY(state.range(0) - 1, state.range(1) - 1));
    ChunkSpanXY<double> src(src_data.data(), dom);
    ChunkSpanXY<double> dst(dst_data.data(), dom);
    for (auto _ : state) {
        ddc::deepcopy(dst, src);
    }
    state.SetBytesProcessed(
            int64_t(state.iterations())
            * int64_t(state.range(0) * state.range(1) * sizeof(double)));
}

static void deepcopy_subchunk_2d(benchmark::State& state)
{
    std::vector<double> src_data(state.range(0) * state.range(1), 0.0);
    std::vector<double> dst_data(state.range(0) * state.range(1), -1.0);
    DDomXY const dom(DElemXY(0, 0), DVectXY(state.range(0) - 1, state.range(1) - 1));
    ChunkSpanXY<double> src(src_data.data(), dom);
    ChunkSpanXY<double> dst(dst_data.data(), dom);
    for (auto _ : state) {
        for (DElemX i : ddc::select<DDimX>(dom)) {
            auto&& dst_i = dst[i];
            auto&& src_i = src[i];
            ddc::deepcopy(dst_i, src_i);
        }
    }
    state.SetBytesProcessed(
            int64_t(state.iterations())
            * int64_t(state.range(0) * state.range(1) * sizeof(double)));
}


// 1D
BENCHMARK(memcpy_1d)->Arg(small_dim1_1D);
BENCHMARK(deepcopy_1d)->Arg(small_dim1_1D);
BENCHMARK(memcpy_1d)->Arg(large_dim1_1D);
BENCHMARK(deepcopy_1d)->Arg(large_dim1_1D);

// 2D
BENCHMARK(memcpy_2d)->Args({small_dim1_2D, small_dim2_2D});
BENCHMARK(deepcopy_subchunk_2d)->Args({small_dim1_2D, small_dim2_2D});
BENCHMARK(memcpy_2d)->Args({large_dim1_2D, large_dim2_2D});
BENCHMARK(deepcopy_2d)->Args({large_dim1_2D, large_dim2_2D});
BENCHMARK(deepcopy_subchunk_2d)->Args({large_dim1_2D, large_dim2_2D});

int main(int argc, char** argv)
{
    ::benchmark::Initialize(&argc, argv);
    if (::benchmark::ReportUnrecognizedArguments(argc, argv)) {
        return 1;
    }
    {
        ddc::ScopeGuard const guard;
        ::benchmark::RunSpecifiedBenchmarks();
    }
    ::benchmark::Shutdown();
    return 0;
}
