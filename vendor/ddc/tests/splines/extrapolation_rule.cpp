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
#include "evaluator_2d.hpp"
#include "polynomial_evaluator.hpp"
#include "spline_error_bounds.hpp"

#if defined(BC_PERIODIC)
struct DimX
{
    static constexpr bool PERIODIC = false;
};

struct DimY
{
    static constexpr bool PERIODIC = true;
};
#elif defined(BC_GREVILLE)

struct DimX
{
    static constexpr bool PERIODIC = false;
};

struct DimY
{
    static constexpr bool PERIODIC = false;
};
#endif

static constexpr std::size_t s_degree = DEGREE;

#if defined(BC_PERIODIC)
static constexpr ddc::BoundCond s_bcl1 = ddc::BoundCond::GREVILLE;
static constexpr ddc::BoundCond s_bcr1 = ddc::BoundCond::GREVILLE;
static constexpr ddc::BoundCond s_bcl2 = ddc::BoundCond::PERIODIC;
static constexpr ddc::BoundCond s_bcr2 = ddc::BoundCond::PERIODIC;
#elif defined(BC_GREVILLE)
static constexpr ddc::BoundCond s_bcl1 = ddc::BoundCond::GREVILLE;
static constexpr ddc::BoundCond s_bcr1 = ddc::BoundCond::GREVILLE;
static constexpr ddc::BoundCond s_bcl2 = ddc::BoundCond::GREVILLE;
static constexpr ddc::BoundCond s_bcr2 = ddc::BoundCond::GREVILLE;
#endif

template <typename BSpX>
using GrevillePoints = ddc::
        GrevilleInterpolationPoints<BSpX, ddc::BoundCond::GREVILLE, ddc::BoundCond::GREVILLE>;

#if defined(BSPLINES_TYPE_UNIFORM)
template <typename X>
using BSplines = ddc::UniformBSplines<X, s_degree>;

// Gives discrete dimension. In the dimension of interest, it is deduced from the BSplines type. In the other dimensions, it has to be newly defined. In practice both types coincide in the test, but it may not be the case.
template <typename X, typename I1, typename I2>
using IDim = std::conditional_t<
        std::is_same_v<X, I1> || std::is_same_v<X, I2>,
        typename GrevillePoints<BSplines<X>>::interpolation_mesh_type,
        ddc::UniformPointSampling<X>>;

#elif defined(BSPLINES_TYPE_NON_UNIFORM)
template <typename X>
using BSplines = ddc::NonUniformBSplines<X, s_degree>;

template <typename X, typename I1, typename I2>
using IDim = std::conditional_t<
        std::is_same_v<X, I1> || std::is_same_v<X, I2>,
        typename GrevillePoints<BSplines<X>>::interpolation_mesh_type,
        ddc::NonUniformPointSampling<X>>;
#endif

#if defined(BC_PERIODIC)
template <typename IDim1, typename IDim2>
using evaluator_type = Evaluator2D::
        Evaluator<CosineEvaluator::Evaluator<IDim1>, CosineEvaluator::Evaluator<IDim2>>;
#else
template <typename IDim1, typename IDim2>
using evaluator_type = Evaluator2D::Evaluator<
        PolynomialEvaluator::Evaluator<IDim1, s_degree>,
        CosineEvaluator::Evaluator<IDim2>>;
#endif

template <typename... IDimX>
using Index = ddc::DiscreteElement<IDimX...>;
template <typename... IDimX>
using DVect = ddc::DiscreteVector<IDimX...>;
template <typename... X>
using Coord = ddc::Coordinate<X...>;

// Extract batch dimensions from IDim (remove dimension of interest). Usefull
template <typename I1, typename I2, typename... X>
using BatchDims = ddc::type_seq_remove_t<ddc::detail::TypeSeq<X...>, ddc::detail::TypeSeq<I1, I2>>;

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
static constexpr double dx(std::size_t ncells)
{
    return (xN<X>() - x0<X>()) / ncells;
}

// Templated function giving break points of mesh in given dimension for non-uniform case.
template <typename X>
static std::vector<Coord<X>> breaks(std::size_t ncells)
{
    std::vector<Coord<X>> out(ncells + 1);
    for (std::size_t i(0); i < ncells + 1; ++i) {
        out[i] = x0<X>() + i * dx<X>(ncells);
    }
    return out;
}

// Helper to initialize space
template <class IDimI1, class IDimI2, class T>
struct DimsInitializer;

template <class IDimI1, class IDimI2, class... IDimX>
struct DimsInitializer<IDimI1, IDimI2, ddc::detail::TypeSeq<IDimX...>>
{
    void operator()(std::size_t const ncells)
    {
#if defined(BSPLINES_TYPE_UNIFORM)
        (ddc::init_discrete_space(IDimX::
                                          init(x0<typename IDimX::continuous_dimension_type>(),
                                               xN<typename IDimX::continuous_dimension_type>(),
                                               DVect<IDimX>(ncells))),
         ...);
        ddc::init_discrete_space<BSplines<typename IDimI1::continuous_dimension_type>>(
                x0<typename IDimI1::continuous_dimension_type>(),
                xN<typename IDimI1::continuous_dimension_type>(),
                ncells);
        ddc::init_discrete_space<BSplines<typename IDimI2::continuous_dimension_type>>(
                x0<typename IDimI2::continuous_dimension_type>(),
                xN<typename IDimI2::continuous_dimension_type>(),
                ncells);
#elif defined(BSPLINES_TYPE_NON_UNIFORM)
        (ddc::init_discrete_space<IDimX>(breaks<typename IDimX::continuous_dimension_type>(ncells)),
         ...);
        ddc::init_discrete_space<BSplines<typename IDimI1::continuous_dimension_type>>(
                breaks<typename IDimI1::continuous_dimension_type>(ncells));
        ddc::init_discrete_space<BSplines<typename IDimI2::continuous_dimension_type>>(
                breaks<typename IDimI2::continuous_dimension_type>(ncells));
#endif
        ddc::init_discrete_space<IDimI1>(
                GrevillePoints<
                        BSplines<typename IDimI1::continuous_dimension_type>>::get_sampling());
        ddc::init_discrete_space<IDimI2>(
                GrevillePoints<
                        BSplines<typename IDimI2::continuous_dimension_type>>::get_sampling());
    }
};

// Checks that when evaluating the spline at interpolation points one
// recovers values that were used to build the spline
template <typename ExecSpace, typename MemorySpace, typename I1, typename I2, typename... X>
static void ExtrapolationRuleSplineTest()
{
    // Instantiate execution spaces and initialize spaces
    Kokkos::DefaultHostExecutionSpace const host_exec_space;
    ExecSpace const exec_space;
    std::size_t constexpr ncells = 10;
    DimsInitializer<
            IDim<I1, I1, I2>,
            IDim<I2, I1, I2>,
            BatchDims<IDim<I1, I1, I2>, IDim<I2, I1, I2>, IDim<X, I1, I2>...>>
            dims_initializer;
    dims_initializer(ncells);

    // Create the values domain (mesh)
    auto interpolation_domain = ddc::DiscreteDomain<IDim<I1, I1, I2>, IDim<I2, I1, I2>>(
            GrevillePoints<BSplines<I1>>::get_domain(),
            GrevillePoints<BSplines<I2>>::get_domain());
    auto const dom_vals_tmp = ddc::DiscreteDomain<IDim<X, void, void>...>(
            ddc::DiscreteDomain<IDim<
                    X,
                    void,
                    void>>(Index<IDim<X, void, void>>(0), DVect<IDim<X, void, void>>(ncells))...);
    ddc::DiscreteDomain<IDim<X, I1, I2>...> const dom_vals
            = ddc::replace_dim_of<IDim<I1, void, void>, IDim<I1, I1, I2>>(
                    ddc::replace_dim_of<
                            IDim<I2, void, void>,
                            IDim<I2, I1, I2>>(dom_vals_tmp, interpolation_domain),
                    interpolation_domain);

    // Create a SplineBuilder over BSplines<I> and batched along other dimensions using some boundary conditions
    ddc::SplineBuilder2D<
            ExecSpace,
            MemorySpace,
            BSplines<I1>,
            BSplines<I2>,
            IDim<I1, I1, I2>,
            IDim<I2, I1, I2>,
            s_bcl1,
            s_bcr1,
            s_bcl2,
            s_bcr2,
            ddc::SplineSolver::GINKGO,
            IDim<X, I1, I2>...>
            spline_builder(dom_vals);

    // Compute usefull domains (dom_interpolation, dom_batch, dom_bsplines and dom_spline)
    ddc::DiscreteDomain<IDim<I1, I1, I2>, IDim<I2, I1, I2>> const dom_interpolation
            = spline_builder.interpolation_domain();
    auto const dom_spline = spline_builder.spline_domain();

    // Allocate and fill a chunk containing values to be passed as input to spline_builder. Those are values of cosine along interest dimension duplicated along batch dimensions
    ddc::Chunk vals1_cpu_alloc(
            dom_interpolation,
            ddc::KokkosAllocator<double, Kokkos::DefaultHostExecutionSpace::memory_space>());
    ddc::ChunkSpan vals1_cpu = vals1_cpu_alloc.span_view();
    evaluator_type<IDim<I1, I1, I2>, IDim<I2, I1, I2>> evaluator(dom_interpolation);
    evaluator(vals1_cpu);
    ddc::Chunk vals1_alloc(dom_interpolation, ddc::KokkosAllocator<double, MemorySpace>());
    ddc::ChunkSpan vals1 = vals1_alloc.span_view();
    ddc::deepcopy(vals1, vals1_cpu);

    ddc::Chunk vals_alloc(dom_vals, ddc::KokkosAllocator<double, MemorySpace>());
    ddc::ChunkSpan vals = vals_alloc.span_view();
    ddc::for_each(
            ddc::policies::policy(exec_space),
            vals.domain(),
            KOKKOS_LAMBDA(Index<IDim<X, I1, I2>...> const e) {
                vals(e) = vals1(ddc::select<IDim<I1, I1, I2>, IDim<I2, I1, I2>>(e));
            });

    // Instantiate chunk of spline coefs to receive output of spline_builder
    ddc::Chunk coef_alloc(dom_spline, ddc::KokkosAllocator<double, MemorySpace>());
    ddc::ChunkSpan coef = coef_alloc.span_view();

    // Finally compute the spline by filling `coef`
    spline_builder(coef, vals.span_cview());

    // Instantiate a SplineEvaluator over interest dimension and batched along other dimensions
    ddc::SplineEvaluator2D<
            ExecSpace,
            MemorySpace,
            BSplines<I1>,
            BSplines<I2>,
            IDim<I1, I1, I2>,
            IDim<I2, I1, I2>,
#if defined(ER_NULL)
            ddc::NullExtrapolationRule,
            ddc::NullExtrapolationRule,
#elif defined(ER_CONSTANT)
            ddc::ConstantExtrapolationRule<I1, I2>,
            ddc::ConstantExtrapolationRule<I1, I2>,
#endif
#if defined(BC_PERIODIC)
            ddc::PeriodicExtrapolationRule<I2>,
            ddc::PeriodicExtrapolationRule<I2>,
#else
#if defined(ER_NULL)
            ddc::NullExtrapolationRule,
            ddc::NullExtrapolationRule,
#elif defined(ER_CONSTANT)
            ddc::ConstantExtrapolationRule<I2, I1>,
            ddc::ConstantExtrapolationRule<I2, I1>,
#endif
#endif

            IDim<X, I1, I2>...>
            spline_evaluator_batched(
                    coef.domain(),
#if defined(ER_NULL)
                    ddc::NullExtrapolationRule(),
                    ddc::NullExtrapolationRule(),
#if defined(BC_PERIODIC)
                    ddc::PeriodicExtrapolationRule<I2>(),
                    ddc::PeriodicExtrapolationRule<I2>()
#else
                    ddc::NullExtrapolationRule(),
                    ddc::NullExtrapolationRule()
#endif
#elif defined(ER_CONSTANT)
#if defined(BC_PERIODIC)
                    ddc::ConstantExtrapolationRule<I1, I2>(x0<I1>()),
                    ddc::ConstantExtrapolationRule<I1, I2>(xN<I1>()),
                    ddc::PeriodicExtrapolationRule<I2>(),
                    ddc::PeriodicExtrapolationRule<I2>()
#else
                    ddc::ConstantExtrapolationRule<I1, I2>(x0<I1>(), x0<I2>(), xN<I2>()),
                    ddc::ConstantExtrapolationRule<I1, I2>(xN<I1>(), x0<I2>(), xN<I2>()),
                    ddc::ConstantExtrapolationRule<I2, I1>(x0<I2>(), x0<I1>(), xN<I1>()),
                    ddc::ConstantExtrapolationRule<I2, I1>(xN<I2>(), x0<I1>(), xN<I1>())
#endif
#endif
            );

    // Instantiate chunk of coordinates of dom_interpolation
    ddc::Chunk coords_eval_alloc(dom_vals, ddc::KokkosAllocator<Coord<X...>, MemorySpace>());
    ddc::ChunkSpan coords_eval = coords_eval_alloc.span_view();
    ddc::for_each(
            ddc::policies::policy(exec_space),
            coords_eval.domain(),
            KOKKOS_LAMBDA(Index<IDim<X, I1, I2>...> const e) {
                coords_eval(e) = ddc::coordinate(e);
                // Set coords_eval outside of the domain
                ddc::get<I1>(coords_eval(e))
                        = xN<I1>() + (ddc::select<I1>(ddc::coordinate(e)) - x0<I1>()) + 1;
#if defined(BC_GREVILLE)
                ddc::get<I2>(coords_eval(e))
                        = xN<I2>() + (ddc::select<I2>(ddc::coordinate(e)) - x0<I2>()) + 1;
#endif
            });


    // Instantiate chunks to receive outputs of spline_evaluator
    ddc::Chunk spline_eval_alloc(dom_vals, ddc::KokkosAllocator<double, MemorySpace>());
    ddc::ChunkSpan spline_eval = spline_eval_alloc.span_view();

    // Call spline_evaluator on the same mesh we started with
    spline_evaluator_batched(spline_eval, coords_eval.span_cview(), coef.span_cview());

    // Checking errors (we recover the initial values)
    double max_norm_error = ddc::transform_reduce(
            ddc::policies::policy(exec_space),
            spline_eval.domain(),
            0.,
            ddc::reducer::max<double>(),
            KOKKOS_LAMBDA(Index<IDim<X, I1, I2>...> const e) {
#if defined(ER_NULL)
                return Kokkos::abs(spline_eval(e));
#elif defined(ER_CONSTANT)
                typename decltype(ddc::remove_dims_of(
                        vals.domain(),
                        vals.template domain<IDim<I1, I1, I2>>()))::discrete_element_type
                        e_without_interest(e);
                typename decltype(ddc::remove_dims_of(
                        vals.domain(),
                        vals.template domain<IDim<I1, I1, I2>, IDim<I2, I1, I2>>()))::
                        discrete_element_type e_batch(e);
                double tmp;
                if (ddc::select<I2>(coords_eval(e)) > xN<I2>()) {
                    tmp = vals(ddc::DiscreteElement<IDim<X, I1, I2>...>(
                            vals.template domain<IDim<I1, I1, I2>>().back(),
                            vals.template domain<IDim<I2, I1, I2>>().back(),
                            e_batch));
                } else {
                    tmp = vals(ddc::DiscreteElement<IDim<X, I1, I2>...>(
                            vals.template domain<IDim<I1, I1, I2>>().back(),
                            e_without_interest));
                }
                return Kokkos::abs(spline_eval(e) - tmp);
#endif
            });

    double const max_norm = evaluator.max_norm();

    EXPECT_LE(max_norm_error, 1.0e-14 * max_norm);
}

#if defined(ER_NULL) && defined(BC_PERIODIC) && defined(BSPLINES_TYPE_UNIFORM)
#define SUFFIX(name) name##Null##Periodic##Uniform
#elif defined(ER_NULL) && defined(BC_PERIODIC) && defined(BSPLINES_TYPE_NON_UNIFORM)
#define SUFFIX(name) name##Null##Periodic##NonUniform
#elif defined(ER_NULL) && defined(BC_GREVILLE) && defined(BSPLINES_TYPE_UNIFORM)
#define SUFFIX(name) name##Null##Greville##Uniform
#elif defined(ER_NULL) && defined(BC_GREVILLE) && defined(BSPLINES_TYPE_NON_UNIFORM)
#define SUFFIX(name) name##Null##Greville##NonUniform
#elif defined(ER_CONSTANT) && defined(BC_PERIODIC) && defined(BSPLINES_TYPE_UNIFORM)
#define SUFFIX(name) name##Constant##Periodic##Uniform
#elif defined(ER_CONSTANT) && defined(BC_PERIODIC) && defined(BSPLINES_TYPE_NON_UNIFORM)
#define SUFFIX(name) name##Constant##Periodic##NonUniform
#elif defined(ER_CONSTANT) && defined(BC_GREVILLE) && defined(BSPLINES_TYPE_UNIFORM)
#define SUFFIX(name) name##Constant##Greville##Uniform
#elif defined(ER_CONSTANT) && defined(BC_GREVILLE) && defined(BSPLINES_TYPE_NON_UNIFORM)
#define SUFFIX(name) name##Constant##Greville##NonUniform
#endif

TEST(SUFFIX(ExtrapolationRuleSplineHost), 2DXY)
{
    ExtrapolationRuleSplineTest<
            Kokkos::DefaultHostExecutionSpace,
            Kokkos::DefaultHostExecutionSpace::memory_space,
            DimX,
            DimY,
            DimX,
            DimY>();
}

TEST(SUFFIX(ExtrapolationRuleSplineDevice), 2DXY)
{
    ExtrapolationRuleSplineTest<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            DimX,
            DimY,
            DimX,
            DimY>();
}
