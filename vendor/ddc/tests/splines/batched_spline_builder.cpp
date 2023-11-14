#include <algorithm>
#include <array>
#include <cmath>
#include <iosfwd>
#include <vector>

#include <experimental/mdspan>

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <gtest/gtest.h>

#include "ddc/coordinate.hpp"
#include "ddc/detail/macros.hpp"
#include "ddc/discrete_domain.hpp"
#include "ddc/for_each.hpp"
#include "ddc/uniform_point_sampling.hpp"

#include "cosine_evaluator.hpp"
#include "polynomial_evaluator.hpp"
#include "spline_error_bounds.hpp"

struct DimX
{
    static constexpr bool PERIODIC = true;
};

struct DimY
{
    static constexpr bool PERIODIC = true;
};

struct DimZ
{
    static constexpr bool PERIODIC = true;
};

struct DimT
{
    static constexpr bool PERIODIC = true;
};

static constexpr std::size_t s_degree_x = DEGREE_X;

template <typename BSpX>
using GrevillePoints = ddc::
        GrevilleInterpolationPoints<BSpX, ddc::BoundCond::PERIODIC, ddc::BoundCond::PERIODIC>;

#if defined(BSPLINES_TYPE_UNIFORM)
template <typename X>
using BSplines = ddc::UniformBSplines<X, s_degree_x>;

// Gives discrete dimension. In the dimension of interest, it is deduced from the BSplines type. In the other dimensions, it has to be newly defined. In practice both types coincide in the test, but it may not be the case.
template <typename X, typename I>
using IDim = std::conditional_t<
        std::is_same_v<X, I>,
        typename GrevillePoints<BSplines<X>>::interpolation_mesh_type,
        ddc::UniformPointSampling<X>>;

#elif defined(BSPLINES_TYPE_NON_UNIFORM)
template <typename X>
using BSplines = ddc::NonUniformBSplines<X, s_degree_x>;

template <typename X, typename I>
using IDim = std::conditional_t<
        std::is_same_v<X, I>,
        typename GrevillePoints<BSplines<X>>::interpolation_mesh_type,
        ddc::NonUniformPointSampling<X>>;
#endif
template <typename IDimX>
using evaluator_type = CosineEvaluator::Evaluator<IDimX>;

template <typename... IDimX>
using Index = ddc::DiscreteElement<IDimX...>;
template <typename... IDimX>
using DVect = ddc::DiscreteVector<IDimX...>;
template <typename... X>
using Coord = ddc::Coordinate<X...>;

// Extract batch dimensions from IDim (remove dimension of interest). Usefull
template <typename I, typename... X>
using BatchDims = ddc::type_seq_remove_t<ddc::detail::TypeSeq<X...>, ddc::detail::TypeSeq<I>>;

// Templated function giving first coordinate of the mesh in given dimension.
template <typename X>
static constexpr Coord<X> x0()
{
    return Coord<X>(0.);
}

// Templated function giving last coordinate of the mesh in given dimension.
template <typename X>
static constexpr Coord<X> xN()
{
    return Coord<X>(1.);
}

// Templated function giving step of the mesh in given dimension.
template <typename X>
static constexpr double dx(double ncells)
{
    return (xN<X>() - x0<X>()) / ncells;
}

// Templated function giving break points of mesh in given dimension for non-uniform case.
template <typename X>
static constexpr std::vector<Coord<X>> breaks(double ncells)
{
    std::vector<Coord<X>> out(ncells + 1);
    for (int i(0); i < ncells + 1; ++i) {
        out[i] = x0<X>() + i * dx<X>(ncells);
    }
    return out;
}

// Helper to initialize space
template <class IDimI, class T>
struct DimsInitializer;

template <class IDimI, class... IDimX>
struct DimsInitializer<IDimI, ddc::detail::TypeSeq<IDimX...>>
{
    void operator()(std::size_t const ncells)
    {
#if defined(BSPLINES_TYPE_UNIFORM)
        (ddc::init_discrete_space(IDimX::
                                          init(x0<typename IDimX::continuous_dimension_type>(),
                                               xN<typename IDimX::continuous_dimension_type>(),
                                               DVect<IDimX>(ncells))),
         ...);
        ddc::init_discrete_space<BSplines<typename IDimI::continuous_dimension_type>>(
                x0<typename IDimI::continuous_dimension_type>(),
                xN<typename IDimI::continuous_dimension_type>(),
                ncells);
#elif defined(BSPLINES_TYPE_NON_UNIFORM)
        (ddc::init_discrete_space<IDimX>(breaks<typename IDimX::continuous_dimension_type>(ncells)),
         ...);
        ddc::init_discrete_space<BSplines<typename IDimI::continuous_dimension_type>>(
                breaks<typename IDimI::continuous_dimension_type>(ncells));
#endif
        ddc::init_discrete_space<IDimI>(
                GrevillePoints<
                        BSplines<typename IDimI::continuous_dimension_type>>::get_sampling());
    }
};

// Checks that when evaluating the spline at interpolation points one
// recovers values that were used to build the spline
template <typename ExecSpace, typename MemorySpace, typename I, typename... X>
static void BatchedPeriodicSplineTest()
{
    // Instantiate execution spaces and initialize spaces
    Kokkos::DefaultHostExecutionSpace host_exec_space = Kokkos::DefaultHostExecutionSpace();
    ExecSpace exec_space = ExecSpace();

    std::size_t constexpr ncells = 10;
    DimsInitializer<IDim<I, I>, BatchDims<IDim<I, I>, IDim<X, I>...>> dims_initializer;
    dims_initializer(ncells);

    // Create the values domain (mesh)
    ddc::DiscreteDomain<IDim<X, I>...> const dom_vals = ddc::DiscreteDomain<IDim<X, I>...>(
            (std::is_same_v<X, I>
                     ? GrevillePoints<BSplines<X>>::get_domain()
                     : ddc::DiscreteDomain<
                             IDim<X, I>>(Index<IDim<X, I>>(0), DVect<IDim<X, I>>(ncells)))...);

    // Create a SplineBuilderBatched over BSplines<I> and batched along other dimensions using some boundary conditions
    ddc::SplineBuilderBatched<
            ddc::SplineBuilder<
                    ExecSpace,
                    MemorySpace,
                    BSplines<I>,
                    IDim<I, I>,
                    ddc::BoundCond::PERIODIC,
                    ddc::BoundCond::PERIODIC>,
            IDim<X, I>...>
            spline_builder(dom_vals);

    // Compute usefull domains (dom_interpolation, dom_batch, dom_bsplines and dom_spline)
    ddc::DiscreteDomain<IDim<I, I>> const dom_interpolation = spline_builder.interpolation_domain();
    auto const dom_batch = spline_builder.batch_domain();
    ddc::DiscreteDomain<BSplines<I>> const dom_bsplines = spline_builder.bsplines_domain();
    auto const dom_spline = spline_builder.spline_domain();

    // Allocate and fill a chunk containing values to be passed as input to spline_builder. Those are values of cosine along interest dimension duplicated along batch dimensions
    ddc::Chunk vals1_cpu_alloc(
            dom_interpolation,
            ddc::KokkosAllocator<double, Kokkos::DefaultHostExecutionSpace::memory_space>());
    ddc::ChunkSpan vals1_cpu = vals1_cpu_alloc.span_view();
    evaluator_type<IDim<I, I>> evaluator(dom_interpolation);
    evaluator(vals1_cpu);
    ddc::Chunk vals1_alloc(dom_interpolation, ddc::KokkosAllocator<double, MemorySpace>());
    ddc::ChunkSpan vals1 = vals1_alloc.span_view();
    ddc::deepcopy(vals1, vals1_cpu);

    ddc::Chunk vals_alloc(dom_vals, ddc::KokkosAllocator<double, MemorySpace>());
    ddc::ChunkSpan vals = vals_alloc.span_view();
    ddc::for_each(
            ddc::policies::policy(exec_space),
            vals.domain(),
            DDC_LAMBDA(Index<IDim<X, I>...> const e) {
                vals(e) = vals1(ddc::select<IDim<I, I>>(e));
            });

    // Instantiate chunk of spline coefs to receive output of spline_builder
    ddc::Chunk coef_alloc(dom_spline, ddc::KokkosAllocator<double, MemorySpace>());
    ddc::ChunkSpan coef = coef_alloc.span_view();

    // Finally compute the spline by filling `coef`
    spline_builder(coef, vals);

    // Instantiate a SplineEvaluator over interest dimension and batched along other dimensions
    ddc::SplineEvaluatorBatched<
            ddc::SplineEvaluator<ExecSpace, MemorySpace, BSplines<I>, IDim<I, I>>,
            IDim<X, I>...>
            spline_evaluator_batched(
                    coef.domain(),
                    ddc::g_null_boundary<BSplines<I>>,
                    ddc::g_null_boundary<BSplines<I>>);

    // Instantiate chunk of coordinates of dom_interpolation
    ddc::Chunk coords_eval_alloc(dom_vals, ddc::KokkosAllocator<Coord<X...>, MemorySpace>());
    ddc::ChunkSpan coords_eval = coords_eval_alloc.span_view();
    ddc::for_each(
            ddc::policies::policy(exec_space),
            coords_eval.domain(),
            DDC_LAMBDA(Index<IDim<X, I>...> const e) { coords_eval(e) = ddc::coordinate(e); });


    // Instantiate chunks to receive outputs of spline_evaluator
    ddc::Chunk spline_eval_alloc(dom_vals, ddc::KokkosAllocator<double, MemorySpace>());
    ddc::ChunkSpan spline_eval = spline_eval_alloc.span_view();
    ddc::Chunk spline_eval_deriv_alloc(dom_vals, ddc::KokkosAllocator<double, MemorySpace>());
    ddc::ChunkSpan spline_eval_deriv = spline_eval_deriv_alloc.span_view();
    ddc::Chunk spline_eval_integrals_alloc(dom_batch, ddc::KokkosAllocator<double, MemorySpace>());
    ddc::ChunkSpan spline_eval_integrals = spline_eval_integrals_alloc.span_view();

    // Call spline_evaluator on the same mesh we started with
    spline_evaluator_batched(spline_eval, coords_eval.span_cview(), coef.span_cview());
    spline_evaluator_batched.deriv(spline_eval_deriv, coords_eval.span_cview(), coef.span_cview());
    spline_evaluator_batched.integrate(spline_eval_integrals, coef.span_cview());

    // Checking errors (we recover the initial values)
    double max_norm_error = ddc::transform_reduce(
            ddc::policies::policy(exec_space),
            spline_eval.domain(),
            0.,
            ddc::reducer::max<double>(),
            DDC_LAMBDA(Index<IDim<X, I>...> const e) {
                return Kokkos::abs(spline_eval(e) - vals(e));
            });

    double max_norm_error_diff = ddc::transform_reduce(
            ddc::policies::policy(exec_space),
            spline_eval_deriv.domain(),
            0.,
            ddc::reducer::max<double>(),
            DDC_LAMBDA(Index<IDim<X, I>...> const e) {
                Coord<I> const x = ddc::coordinate(ddc::select<IDim<I, I>>(e));
                return Kokkos::abs(spline_eval_deriv(e) - evaluator.deriv(x, 1));
            });
    double max_norm_error_integ = ddc::transform_reduce(
            ddc::policies::policy(exec_space),
            spline_eval_integrals.domain(),
            0.,
            ddc::reducer::max<double>(),
            DDC_LAMBDA(typename decltype(spline_builder)::batch_domain_type::
                               discrete_element_type const e) {
                return Kokkos::abs(
                        spline_eval_integrals(e) - evaluator.deriv(xN<I>(), -1)
                        + evaluator.deriv(x0<I>(), -1));
            });

    double const max_norm = evaluator.max_norm();
    double const max_norm_diff = evaluator.max_norm(1);
    double const max_norm_int = evaluator.max_norm(-1);

    SplineErrorBounds<evaluator_type<IDim<I, I>>> error_bounds(evaluator);
    EXPECT_LE(
            max_norm_error,
            std::max(error_bounds.error_bound(dx<I>(ncells), s_degree_x), 1.0e-14 * max_norm));
    EXPECT_LE(
            max_norm_error_diff,
            std::
                    max(error_bounds.error_bound_on_deriv(dx<I>(ncells), s_degree_x),
                        1e-12 * max_norm_diff));
    EXPECT_LE(
            max_norm_error_integ,
            std::
                    max(error_bounds.error_bound_on_int(dx<I>(ncells), s_degree_x),
                        1.0e-14 * max_norm_int));
}

TEST(BatchedPeriodicSplineHost, 1DX)
{
    BatchedPeriodicSplineTest<
            Kokkos::DefaultHostExecutionSpace,
            Kokkos::DefaultHostExecutionSpace::memory_space,
            DimX,
            DimX>();
}

TEST(BatchedPeriodicSplineDevice, 1DX)
{
    BatchedPeriodicSplineTest<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            DimX,
            DimX>();
}

TEST(BatchedPeriodicSplineHost, 2DX)
{
    BatchedPeriodicSplineTest<
            Kokkos::DefaultHostExecutionSpace,
            Kokkos::DefaultHostExecutionSpace::memory_space,
            DimX,
            DimX,
            DimY>();
}

TEST(BatchedPeriodicSplineHost, 2DY)
{
    BatchedPeriodicSplineTest<
            Kokkos::DefaultHostExecutionSpace,
            Kokkos::DefaultHostExecutionSpace::memory_space,
            DimY,
            DimX,
            DimY>();
}

TEST(BatchedPeriodicSplineDevice, 2DX)
{
    BatchedPeriodicSplineTest<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            DimX,
            DimX,
            DimY>();
}

TEST(BatchedPeriodicSplineDevice, 2DY)
{
    BatchedPeriodicSplineTest<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            DimY,
            DimX,
            DimY>();
}

TEST(BatchedPeriodicSplineHost, 3DX)
{
    BatchedPeriodicSplineTest<
            Kokkos::DefaultHostExecutionSpace,
            Kokkos::DefaultHostExecutionSpace::memory_space,
            DimX,
            DimX,
            DimY,
            DimZ>();
}

TEST(BatchedPeriodicSplineHost, 3DY)
{
    BatchedPeriodicSplineTest<
            Kokkos::DefaultHostExecutionSpace,
            Kokkos::DefaultHostExecutionSpace::memory_space,
            DimY,
            DimX,
            DimY,
            DimZ>();
}

TEST(BatchedPeriodicSplineHost, 3DZ)
{
    BatchedPeriodicSplineTest<
            Kokkos::DefaultHostExecutionSpace,
            Kokkos::DefaultHostExecutionSpace::memory_space,
            DimZ,
            DimX,
            DimY,
            DimZ>();
}

TEST(BatchedPeriodicSplineDevice, 3DX)
{
    BatchedPeriodicSplineTest<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            DimX,
            DimX,
            DimY,
            DimZ>();
}

TEST(BatchedPeriodicSplineDevice, 3DY)
{
    BatchedPeriodicSplineTest<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            DimY,
            DimX,
            DimY,
            DimZ>();
}

TEST(BatchedPeriodicSplineDevice, 3DZ)
{
    BatchedPeriodicSplineTest<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            DimZ,
            DimX,
            DimY,
            DimZ>();
}


TEST(BatchedPeriodicSplineHost, 4DX)
{
    BatchedPeriodicSplineTest<
            Kokkos::DefaultHostExecutionSpace,
            Kokkos::DefaultHostExecutionSpace::memory_space,
            DimX,
            DimX,
            DimY,
            DimZ,
            DimT>();
}

TEST(BatchedPeriodicSplineHost, 4DY)
{
    BatchedPeriodicSplineTest<
            Kokkos::DefaultHostExecutionSpace,
            Kokkos::DefaultHostExecutionSpace::memory_space,
            DimY,
            DimX,
            DimY,
            DimZ,
            DimT>();
}

TEST(BatchedPeriodicSplineHost, 4DZ)
{
    BatchedPeriodicSplineTest<
            Kokkos::DefaultHostExecutionSpace,
            Kokkos::DefaultHostExecutionSpace::memory_space,
            DimZ,
            DimX,
            DimY,
            DimZ,
            DimT>();
}

TEST(BatchedPeriodicSplineHost, 4DT)
{
    BatchedPeriodicSplineTest<
            Kokkos::DefaultHostExecutionSpace,
            Kokkos::DefaultHostExecutionSpace::memory_space,
            DimT,
            DimX,
            DimY,
            DimZ,
            DimT>();
}

TEST(BatchedPeriodicSplineDevice, 4DX)
{
    BatchedPeriodicSplineTest<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            DimX,
            DimX,
            DimY,
            DimZ,
            DimT>();
}

TEST(BatchedPeriodicSplineDevice, 4DY)
{
    BatchedPeriodicSplineTest<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            DimY,
            DimX,
            DimY,
            DimZ,
            DimT>();
}

TEST(BatchedPeriodicSplineDevice, 4DZ)
{
    BatchedPeriodicSplineTest<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            DimZ,
            DimX,
            DimY,
            DimZ,
            DimT>();
}

TEST(BatchedPeriodicSplineDevice, 4DT)
{
    BatchedPeriodicSplineTest<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            DimT,
            DimX,
            DimY,
            DimZ,
            DimT>();
}
