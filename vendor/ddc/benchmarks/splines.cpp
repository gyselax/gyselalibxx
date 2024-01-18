// SPDX-License-Identifier: MIT
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <iosfwd>
#include <vector>

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <benchmark/benchmark.h>

namespace {

static constexpr std::size_t s_degree_x = 3;

struct X
{
    static constexpr bool PERIODIC = true;
};

using BSplinesX = ddc::UniformBSplines<X, s_degree_x>;
using GrevillePoints = ddc::
        GrevilleInterpolationPoints<BSplinesX, ddc::BoundCond::PERIODIC, ddc::BoundCond::PERIODIC>;
using DDimX = GrevillePoints::interpolation_mesh_type;

struct Y;
using DDimY = ddc::UniformPointSampling<Y>;


} // namespace

// Function to monitor GPU memory asynchronously
void monitorMemoryAsync(std::mutex& mutex, bool& monitorFlag, size_t& maxUsedMem)
{
    size_t freeMem = 0;
    size_t totalMem = 0;
    while (monitorFlag) {
        std::this_thread::sleep_for(std::chrono::milliseconds(1)); // Adjust the interval as needed

        // Acquire a lock to ensure thread safety when accessing CUDA functions
        std::lock_guard<std::mutex> lock(mutex);

#if defined(__CUDACC__)
        cudaMemGetInfo(&freeMem, &totalMem);
#endif
        maxUsedMem = std::max(maxUsedMem, totalMem - freeMem);
    }
}

static void characteristics_advection(benchmark::State& state)
{
    size_t freeMem = 0;
    size_t totalMem = 0;
#if defined(__CUDACC__)
    cudaMemGetInfo(&freeMem, &totalMem);
#endif
    size_t initUsedMem
            = totalMem
              - freeMem; // cudaMemGetInfo gives GPU total memory occupation, we consider that other users of the GPU have constant occupancy and substract it.
    size_t maxUsedMem = initUsedMem;

    bool monitorFlag = true;
    std::mutex mutex;
    // Create a thread to monitor GPU memory asynchronously
    std::thread monitorThread(
            monitorMemoryAsync,
            std::ref(mutex),
            std::ref(monitorFlag),
            std::ref(maxUsedMem));

    ddc::init_discrete_space<
            BSplinesX>(ddc::Coordinate<X>(-1.), ddc::Coordinate<X>(1.), state.range(0));
    ddc::init_discrete_space<DDimX>(ddc::GrevilleInterpolationPoints<
                                    BSplinesX,
                                    ddc::BoundCond::PERIODIC,
                                    ddc::BoundCond::PERIODIC>::get_sampling());
    ddc::DiscreteDomain<DDimY> y_domain
            = ddc::init_discrete_space(DDimY::
                                               init(ddc::Coordinate<Y>(-1.),
                                                    ddc::Coordinate<Y>(1.),
                                                    ddc::DiscreteVector<DDimY>(state.range(1))));

    auto const x_domain = ddc::GrevilleInterpolationPoints<
            BSplinesX,
            ddc::BoundCond::PERIODIC,
            ddc::BoundCond::PERIODIC>::get_domain();
    ddc::Chunk density_alloc(
            ddc::DiscreteDomain<DDimX, DDimY>(x_domain, y_domain),
            ddc::DeviceAllocator<double>());
    ddc::ChunkSpan const density = density_alloc.span_view();
    // Initialize the density on the main domain
    ddc::DiscreteDomain<DDimX, DDimY> x_mesh
            = ddc::DiscreteDomain<DDimX, DDimY>(x_domain, y_domain);
    ddc::for_each(
            ddc::policies::parallel_device,
            x_mesh,
            KOKKOS_LAMBDA(ddc::DiscreteElement<DDimX, DDimY> const ixy) {
                double const x = ddc::coordinate(ddc::select<DDimX>(ixy));
                double const y = ddc::coordinate(ddc::select<DDimY>(ixy));
                density(ixy) = 9.999 * Kokkos::exp(-(x * x + y * y) / 0.1 / 2);
                // initial_density(ixy) = 9.999 * ((x * x + y * y) < 0.25);
            });
    ddc::SplineBuilderBatched<
            ddc::SplineBuilder<
                    Kokkos::DefaultExecutionSpace,
                    Kokkos::DefaultExecutionSpace::memory_space,
                    BSplinesX,
                    DDimX,
                    ddc::BoundCond::PERIODIC,
                    ddc::BoundCond::PERIODIC>,
            DDimX,
            DDimY>
            spline_builder(x_mesh, state.range(2), state.range(3));
    ddc::SplineEvaluatorBatched<
            ddc::SplineEvaluator<
                    Kokkos::DefaultExecutionSpace,
                    Kokkos::DefaultExecutionSpace::memory_space,
                    BSplinesX,
                    DDimX>,
            DDimX,
            DDimY>
            spline_evaluator(
                    spline_builder.spline_domain(),
                    ddc::g_null_boundary<BSplinesX>,
                    ddc::g_null_boundary<BSplinesX>);
    ddc::Chunk coef_alloc(
            spline_builder.spline_domain(),
            ddc::KokkosAllocator<double, Kokkos::DefaultExecutionSpace::memory_space>());
    ddc::ChunkSpan coef = coef_alloc.span_view();
    ddc::Chunk feet_coords_alloc(
            spline_builder.vals_domain(),
            ddc::KokkosAllocator<
                    ddc::Coordinate<X, Y>,
                    Kokkos::DefaultExecutionSpace::memory_space>());
    ddc::ChunkSpan feet_coords = feet_coords_alloc.span_view();

    for (auto _ : state) {
        Kokkos::Profiling::pushRegion("FeetCharacteristics");
        ddc::for_each(
                ddc::policies::parallel_device,
                feet_coords.domain(),
                KOKKOS_LAMBDA(ddc::DiscreteElement<DDimX, DDimY> const e) {
                    feet_coords(e) = ddc::Coordinate<X, Y>(
                            ddc::coordinate(ddc::select<DDimX>(e))
                                    - ddc::Coordinate<X>(0.0176429863),
                            ddc::coordinate(ddc::select<DDimY>(e)));
                });
        Kokkos::Profiling::popRegion();
        Kokkos::Profiling::pushRegion("SplineBuilder");
        spline_builder(coef, density.span_cview());
        Kokkos::Profiling::popRegion();
        Kokkos::Profiling::pushRegion("SplineEvaluator");
        spline_evaluator(density, feet_coords.span_cview(), coef.span_cview());
        Kokkos::Profiling::popRegion();
    }
    monitorFlag = false;
    monitorThread.join();
    state.SetBytesProcessed(
            int64_t(state.iterations())
            * int64_t(state.range(0) * state.range(1) * sizeof(double)));
    state.counters["gpu_mem_occupancy"] = maxUsedMem - initUsedMem;
    ////////////////////////////////////////////////////
    /// --------------- HUGE WARNING --------------- ///
    /// The following lines are forbidden in a prod- ///
    /// uction code. It is a necessary workaround    ///
    /// which must be used ONLY for Google Benchmark.///
    /// The reason is it acts on underlying global   ///
    /// variables, which is always a bad idea.       ///
    ////////////////////////////////////////////////////
    ddc::detail::g_discrete_space_dual<BSplinesX>.reset();
    ddc::detail::g_discrete_space_dual<BSplinesX::mesh_type>.reset();
    ddc::detail::g_discrete_space_dual<DDimX>.reset();
    ddc::detail::g_discrete_space_dual<DDimY>.reset();
    ////////////////////////////////////////////////////
}

// Tuning : 512 cols and 8 precond on CPU, 16384 cols and 1 precond on GPU

#ifdef KOKKOS_ENABLE_CUDA
std::string chip = "gpu";
int cols_per_chunk_ref = 65535;
unsigned int preconditionner_max_block_size_ref = 1u;
#elif defined(KOKKOS_ENABLE_OPENMP)
std::string chip = "cpu";
int cols_per_chunk_ref = 8192;
unsigned int preconditionner_max_block_size_ref = 32u;
#elif defined(KOKKOS_ENABLE_SERIAL)
std::string chip = "cpu";
int cols_per_chunk_ref = 8192;
unsigned int preconditionner_max_block_size_ref = 32u;
#endif

BENCHMARK(characteristics_advection)
        ->RangeMultiplier(2)
        ->Ranges(
                {{64, 1024},
                 {100, 200000},
                 {cols_per_chunk_ref, cols_per_chunk_ref},
                 {preconditionner_max_block_size_ref, preconditionner_max_block_size_ref}})
        ->MinTime(3)
        ->UseRealTime();
/*
BENCHMARK(characteristics_advection)
        ->RangeMultiplier(2)
        ->Ranges({{64, 1024}, {100000, 100000}, {64,65535}, {preconditionner_max_block_size_ref, preconditionner_max_block_size_ref}})
        ->MinTime(3)->UseRealTime();
*/
/*
BENCHMARK(characteristics_advection)
        ->RangeMultiplier(2)
        ->Ranges({{64, 1024}, {100000, 100000}, {cols_per_chunk_ref, cols_per_chunk_ref}, {1, 32}})
        ->MinTime(3)->UseRealTime();
*/

int main(int argc, char** argv)
{
    ::benchmark::Initialize(&argc, argv);
    ::benchmark::AddCustomContext("chip", chip);
    ::benchmark::AddCustomContext("cols_per_chunk_ref", std::to_string(cols_per_chunk_ref));
    ::benchmark::AddCustomContext(
            "preconditionner_max_block_size_ref",
            std::to_string(preconditionner_max_block_size_ref));
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
