// Copyright (C) The DDC development team, see COPYRIGHT.md file
//
// SPDX-License-Identifier: MIT

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
#else

struct DimX
{
    static constexpr bool PERIODIC = false;
};

struct DimY
{
    static constexpr bool PERIODIC = false;
};

struct DimZ
{
    static constexpr bool PERIODIC = false;
};

struct DimT
{
    static constexpr bool PERIODIC = false;
};
#endif

static constexpr std::size_t s_degree = DEGREE;

#if defined(BC_PERIODIC)
static constexpr ddc::BoundCond s_bcl = ddc::BoundCond::PERIODIC;
static constexpr ddc::BoundCond s_bcr = ddc::BoundCond::PERIODIC;
#elif defined(BC_GREVILLE)
static constexpr ddc::BoundCond s_bcl = ddc::BoundCond::GREVILLE;
static constexpr ddc::BoundCond s_bcr = ddc::BoundCond::GREVILLE;
#elif defined(BC_HERMITE)
static constexpr ddc::BoundCond s_bcl = ddc::BoundCond::HERMITE;
static constexpr ddc::BoundCond s_bcr = ddc::BoundCond::HERMITE;
#endif

template <typename BSpX>
using GrevillePoints = ddc::GrevilleInterpolationPoints<BSpX, s_bcl, s_bcr>;

#if defined(BSPLINES_TYPE_UNIFORM)
template <typename X>
struct BSplines : ddc::UniformBSplines<X, s_degree>
{
};

// Gives discrete dimension. In the dimension of interest, it is deduced from the BSplines type. In the other dimensions, it has to be newly defined. In practice both types coincide in the test, but it may not be the case.
template <class X, bool B>
struct IDim_
    : std::conditional_t<
              B,
              typename GrevillePoints<BSplines<X>>::interpolation_mesh_type,
              ddc::UniformPointSampling<X>>
{
};

template <typename X, typename I1, typename I2>
using IDim = IDim_<X, std::is_same_v<X, I1> || std::is_same_v<X, I2>>;

#elif defined(BSPLINES_TYPE_NON_UNIFORM)
template <typename X>
struct BSplines : ddc::NonUniformBSplines<X, s_degree>
{
};

template <class X>
struct IDim_ : ddc::NonUniformPointSampling<X>
{
};

template <typename X, typename I1, typename I2>
using IDim = IDim_<X>;
#endif

#if defined(BC_PERIODIC)
template <typename IDim1, typename IDim2>
using evaluator_type = Evaluator2D::
        Evaluator<CosineEvaluator::Evaluator<IDim1>, CosineEvaluator::Evaluator<IDim2>>;
#else
template <typename IDim1, typename IDim2>
using evaluator_type = Evaluator2D::Evaluator<
        PolynomialEvaluator::Evaluator<IDim1, s_degree>,
        PolynomialEvaluator::Evaluator<IDim2, s_degree>>;
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
        (ddc::init_discrete_space<IDimX>(IDimX::template init<IDimX>(
                 x0<typename IDimX::continuous_dimension_type>(),
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
                GrevillePoints<BSplines<typename IDimI1::continuous_dimension_type>>::
                        template get_sampling<IDimI1>());
        ddc::init_discrete_space<IDimI2>(
                GrevillePoints<BSplines<typename IDimI2::continuous_dimension_type>>::
                        template get_sampling<IDimI2>());
    }
};

// Checks that when evaluating the spline at interpolation points one
// recovers values that were used to build the spline
template <typename ExecSpace, typename MemorySpace, typename I1, typename I2, typename... X>
static void Batched2dSplineTest()
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
#if defined(BC_HERMITE)
    auto interpolation_domain1 = ddc::DiscreteDomain<IDim<I1, I1, I2>>(
            GrevillePoints<BSplines<I1>>::template get_domain<IDim<I1, I1, I2>>());
    auto interpolation_domain2 = ddc::DiscreteDomain<IDim<I2, I1, I2>>(
            GrevillePoints<BSplines<I2>>::template get_domain<IDim<I2, I1, I2>>());
#endif
    auto interpolation_domain = ddc::DiscreteDomain<IDim<I1, I1, I2>, IDim<I2, I1, I2>>(
            GrevillePoints<BSplines<I1>>::template get_domain<IDim<I1, I1, I2>>(),
            GrevillePoints<BSplines<I2>>::template get_domain<IDim<I2, I1, I2>>());
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

#if defined(BC_HERMITE)
    // Create the derivs domain
    ddc::DiscreteDomain<ddc::Deriv<I1>> const derivs_domain1 = ddc::DiscreteDomain<
            ddc::Deriv<I1>>(Index<ddc::Deriv<I1>>(1), DVect<ddc::Deriv<I1>>(s_degree / 2));
    ddc::DiscreteDomain<ddc::Deriv<I2>> const derivs_domain2 = ddc::DiscreteDomain<
            ddc::Deriv<I2>>(Index<ddc::Deriv<I2>>(1), DVect<ddc::Deriv<I2>>(s_degree / 2));
    ddc::DiscreteDomain<ddc::Deriv<I1>, ddc::Deriv<I2>> const derivs_domain
            = ddc::DiscreteDomain<ddc::Deriv<I1>, ddc::Deriv<I2>>(derivs_domain1, derivs_domain2);

    auto const dom_derivs1
            = ddc::replace_dim_of<IDim<I1, I1, I2>, ddc::Deriv<I1>>(dom_vals, derivs_domain1);
    auto const dom_derivs2
            = ddc::replace_dim_of<IDim<I2, I1, I2>, ddc::Deriv<I2>>(dom_vals, derivs_domain2);
    auto const dom_derivs
            = ddc::replace_dim_of<IDim<I2, I1, I2>, ddc::Deriv<I2>>(dom_derivs1, derivs_domain2);
#endif

    // Create a SplineBuilder over BSplines<I> and batched along other dimensions using some boundary conditions
    ddc::SplineBuilder2D<
            ExecSpace,
            MemorySpace,
            BSplines<I1>,
            BSplines<I2>,
            IDim<I1, I1, I2>,
            IDim<I2, I1, I2>,
            s_bcl,
            s_bcr,
            s_bcl,
            s_bcr,
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
    ddc::parallel_deepcopy(vals1, vals1_cpu);

    ddc::Chunk vals_alloc(dom_vals, ddc::KokkosAllocator<double, MemorySpace>());
    ddc::ChunkSpan vals = vals_alloc.span_view();
    ddc::parallel_for_each(
            exec_space,
            vals.domain(),
            KOKKOS_LAMBDA(Index<IDim<X, I1, I2>...> const e) {
                vals(e) = vals1(ddc::select<IDim<I1, I1, I2>, IDim<I2, I1, I2>>(e));
            });

#if defined(BC_HERMITE)
    // Allocate and fill a chunk containing derivs to be passed as input to spline_builder.
    int constexpr shift = s_degree % 2; // shift = 0 for even order, 1 for odd order
    ddc::Chunk Sderiv1_lhs_alloc(dom_derivs1, ddc::KokkosAllocator<double, MemorySpace>());
    ddc::ChunkSpan Sderiv1_lhs = Sderiv1_lhs_alloc.span_view();
    if (s_bcl == ddc::BoundCond::HERMITE) {
        ddc::Chunk Sderiv1_lhs1_cpu_alloc(
                ddc::DiscreteDomain<
                        ddc::Deriv<I1>,
                        IDim<I2, I1, I2>>(derivs_domain1, interpolation_domain2),
                ddc::HostAllocator<double>());
        ddc::ChunkSpan Sderiv1_lhs1_cpu = Sderiv1_lhs1_cpu_alloc.span_view();
        ddc::for_each(
                Sderiv1_lhs1_cpu.domain(),
                KOKKOS_LAMBDA(ddc::DiscreteElement<ddc::Deriv<I1>, IDim<I2, I1, I2>> const e) {
                    auto deriv_idx = ddc::DiscreteElement<ddc::Deriv<I1>>(e).uid();
                    auto x2 = ddc::coordinate(ddc::DiscreteElement<IDim<I2, I1, I2>>(e));
                    Sderiv1_lhs1_cpu(e) = evaluator.deriv(x0<I1>(), x2, deriv_idx + shift - 1, 0);
                });
        ddc::Chunk Sderiv1_lhs1_alloc(
                ddc::DiscreteDomain<
                        ddc::Deriv<I1>,
                        IDim<I2, I1, I2>>(derivs_domain1, interpolation_domain2),
                ddc::KokkosAllocator<double, MemorySpace>());
        ddc::ChunkSpan Sderiv1_lhs1 = Sderiv1_lhs1_alloc.span_view();
        ddc::parallel_deepcopy(Sderiv1_lhs1, Sderiv1_lhs1_cpu);

        ddc::parallel_for_each(
                exec_space,
                Sderiv1_lhs.domain(),
                KOKKOS_LAMBDA(
                        typename decltype(Sderiv1_lhs.domain())::discrete_element_type const e) {
                    Sderiv1_lhs(e) = Sderiv1_lhs1(ddc::select<ddc::Deriv<I1>, IDim<I2, I1, I2>>(e));
                });
    }

    ddc::Chunk Sderiv1_rhs_alloc(dom_derivs1, ddc::KokkosAllocator<double, MemorySpace>());
    ddc::ChunkSpan Sderiv1_rhs = Sderiv1_rhs_alloc.span_view();
    if (s_bcl == ddc::BoundCond::HERMITE) {
        ddc::Chunk Sderiv1_rhs1_cpu_alloc(
                ddc::DiscreteDomain<
                        ddc::Deriv<I1>,
                        IDim<I2, I1, I2>>(derivs_domain1, interpolation_domain2),
                ddc::HostAllocator<double>());
        ddc::ChunkSpan Sderiv1_rhs1_cpu = Sderiv1_rhs1_cpu_alloc.span_view();
        ddc::for_each(
                Sderiv1_rhs1_cpu.domain(),
                KOKKOS_LAMBDA(ddc::DiscreteElement<ddc::Deriv<I1>, IDim<I2, I1, I2>> const e) {
                    auto deriv_idx = ddc::DiscreteElement<ddc::Deriv<I1>>(e).uid();
                    auto x2 = ddc::coordinate(ddc::DiscreteElement<IDim<I2, I1, I2>>(e));
                    Sderiv1_rhs1_cpu(e) = evaluator.deriv(xN<I1>(), x2, deriv_idx + shift - 1, 0);
                });
        ddc::Chunk Sderiv1_rhs1_alloc(
                ddc::DiscreteDomain<
                        ddc::Deriv<I1>,
                        IDim<I2, I1, I2>>(derivs_domain1, interpolation_domain2),
                ddc::KokkosAllocator<double, MemorySpace>());
        ddc::ChunkSpan Sderiv1_rhs1 = Sderiv1_rhs1_alloc.span_view();
        ddc::parallel_deepcopy(Sderiv1_rhs1, Sderiv1_rhs1_cpu);

        ddc::parallel_for_each(
                exec_space,
                Sderiv1_rhs.domain(),
                KOKKOS_LAMBDA(
                        typename decltype(Sderiv1_rhs.domain())::discrete_element_type const e) {
                    Sderiv1_rhs(e) = Sderiv1_rhs1(ddc::select<ddc::Deriv<I1>, IDim<I2, I1, I2>>(e));
                });
    }

    ddc::Chunk Sderiv2_lhs_alloc(dom_derivs2, ddc::KokkosAllocator<double, MemorySpace>());
    ddc::ChunkSpan Sderiv2_lhs = Sderiv2_lhs_alloc.span_view();
    if (s_bcl == ddc::BoundCond::HERMITE) {
        ddc::Chunk Sderiv2_lhs1_cpu_alloc(
                ddc::DiscreteDomain<
                        IDim<I1, I1, I2>,
                        ddc::Deriv<I2>>(interpolation_domain1, derivs_domain2),
                ddc::HostAllocator<double>());
        ddc::ChunkSpan Sderiv2_lhs1_cpu = Sderiv2_lhs1_cpu_alloc.span_view();
        ddc::for_each(
                Sderiv2_lhs1_cpu.domain(),
                KOKKOS_LAMBDA(ddc::DiscreteElement<IDim<I1, I1, I2>, ddc::Deriv<I2>> const e) {
                    auto x1 = ddc::coordinate(ddc::DiscreteElement<IDim<I1, I1, I2>>(e));
                    auto deriv_idx = ddc::DiscreteElement<ddc::Deriv<I2>>(e).uid();
                    Sderiv2_lhs1_cpu(e) = evaluator.deriv(x1, x0<I2>(), 0, deriv_idx + shift - 1);
                });

        ddc::Chunk Sderiv2_lhs1_alloc(
                ddc::DiscreteDomain<
                        IDim<I1, I1, I2>,
                        ddc::Deriv<I2>>(interpolation_domain1, derivs_domain2),
                ddc::KokkosAllocator<double, MemorySpace>());
        ddc::ChunkSpan Sderiv2_lhs1 = Sderiv2_lhs1_alloc.span_view();
        ddc::parallel_deepcopy(Sderiv2_lhs1, Sderiv2_lhs1_cpu);

        ddc::parallel_for_each(
                exec_space,
                Sderiv2_lhs.domain(),
                KOKKOS_LAMBDA(
                        typename decltype(Sderiv2_lhs.domain())::discrete_element_type const e) {
                    Sderiv2_lhs(e) = Sderiv2_lhs1(ddc::select<IDim<I1, I1, I2>, ddc::Deriv<I2>>(e));
                });
    }

    ddc::Chunk Sderiv2_rhs_alloc(dom_derivs2, ddc::KokkosAllocator<double, MemorySpace>());
    ddc::ChunkSpan Sderiv2_rhs = Sderiv2_rhs_alloc.span_view();
    if (s_bcl == ddc::BoundCond::HERMITE) {
        ddc::Chunk Sderiv2_rhs1_cpu_alloc(
                ddc::DiscreteDomain<
                        IDim<I1, I1, I2>,
                        ddc::Deriv<I2>>(interpolation_domain1, derivs_domain2),
                ddc::HostAllocator<double>());
        ddc::ChunkSpan Sderiv2_rhs1_cpu = Sderiv2_rhs1_cpu_alloc.span_view();
        ddc::for_each(
                Sderiv2_rhs1_cpu.domain(),
                KOKKOS_LAMBDA(ddc::DiscreteElement<IDim<I1, I1, I2>, ddc::Deriv<I2>> const e) {
                    auto x1 = ddc::coordinate(ddc::DiscreteElement<IDim<I1, I1, I2>>(e));
                    auto deriv_idx = ddc::DiscreteElement<ddc::Deriv<I2>>(e).uid();
                    Sderiv2_rhs1_cpu(e) = evaluator.deriv(x1, xN<I2>(), 0, deriv_idx + shift - 1);
                });

        ddc::Chunk Sderiv2_rhs1_alloc(
                ddc::DiscreteDomain<
                        IDim<I1, I1, I2>,
                        ddc::Deriv<I2>>(interpolation_domain1, derivs_domain2),
                ddc::KokkosAllocator<double, MemorySpace>());
        ddc::ChunkSpan Sderiv2_rhs1 = Sderiv2_rhs1_alloc.span_view();
        ddc::parallel_deepcopy(Sderiv2_rhs1, Sderiv2_rhs1_cpu);

        ddc::parallel_for_each(
                exec_space,
                Sderiv2_rhs.domain(),
                KOKKOS_LAMBDA(
                        typename decltype(Sderiv2_rhs.domain())::discrete_element_type const e) {
                    Sderiv2_rhs(e) = Sderiv2_rhs1(ddc::select<IDim<I1, I1, I2>, ddc::Deriv<I2>>(e));
                });
    }

    ddc::Chunk Sderiv_mixed_lhs_lhs_alloc(dom_derivs, ddc::KokkosAllocator<double, MemorySpace>());
    ddc::ChunkSpan Sderiv_mixed_lhs_lhs = Sderiv_mixed_lhs_lhs_alloc.span_view();
    ddc::Chunk Sderiv_mixed_rhs_lhs_alloc(dom_derivs, ddc::KokkosAllocator<double, MemorySpace>());
    ddc::ChunkSpan Sderiv_mixed_rhs_lhs = Sderiv_mixed_rhs_lhs_alloc.span_view();
    ddc::Chunk Sderiv_mixed_lhs_rhs_alloc(dom_derivs, ddc::KokkosAllocator<double, MemorySpace>());
    ddc::ChunkSpan Sderiv_mixed_lhs_rhs = Sderiv_mixed_lhs_rhs_alloc.span_view();
    ddc::Chunk Sderiv_mixed_rhs_rhs_alloc(dom_derivs, ddc::KokkosAllocator<double, MemorySpace>());
    ddc::ChunkSpan Sderiv_mixed_rhs_rhs = Sderiv_mixed_rhs_rhs_alloc.span_view();

    if (s_bcl == ddc::BoundCond::HERMITE) {
        ddc::Chunk Sderiv_mixed_lhs_lhs1_cpu_alloc(derivs_domain, ddc::HostAllocator<double>());
        ddc::ChunkSpan Sderiv_mixed_lhs_lhs1_cpu = Sderiv_mixed_lhs_lhs1_cpu_alloc.span_view();
        ddc::Chunk Sderiv_mixed_rhs_lhs1_cpu_alloc(derivs_domain, ddc::HostAllocator<double>());
        ddc::ChunkSpan Sderiv_mixed_rhs_lhs1_cpu = Sderiv_mixed_rhs_lhs1_cpu_alloc.span_view();
        ddc::Chunk Sderiv_mixed_lhs_rhs1_cpu_alloc(derivs_domain, ddc::HostAllocator<double>());
        ddc::ChunkSpan Sderiv_mixed_lhs_rhs1_cpu = Sderiv_mixed_lhs_rhs1_cpu_alloc.span_view();
        ddc::Chunk Sderiv_mixed_rhs_rhs1_cpu_alloc(derivs_domain, ddc::HostAllocator<double>());
        ddc::ChunkSpan Sderiv_mixed_rhs_rhs1_cpu = Sderiv_mixed_rhs_rhs1_cpu_alloc.span_view();

        for (std::size_t ii = 1;
             ii < (std::size_t)derivs_domain.template extent<ddc::Deriv<I1>>() + 1;
             ++ii) {
            for (std::size_t jj = 1;
                 jj < (std::size_t)derivs_domain.template extent<ddc::Deriv<I2>>() + 1;
                 ++jj) {
                Sderiv_mixed_lhs_lhs1_cpu(
                        typename decltype(derivs_domain)::discrete_element_type(ii, jj))
                        = evaluator.deriv(x0<I1>(), x0<I2>(), ii + shift - 1, jj + shift - 1);
                Sderiv_mixed_rhs_lhs1_cpu(
                        typename decltype(derivs_domain)::discrete_element_type(ii, jj))
                        = evaluator.deriv(xN<I1>(), x0<I2>(), ii + shift - 1, jj + shift - 1);
                Sderiv_mixed_lhs_rhs1_cpu(
                        typename decltype(derivs_domain)::discrete_element_type(ii, jj))
                        = evaluator.deriv(x0<I1>(), xN<I2>(), ii + shift - 1, jj + shift - 1);
                Sderiv_mixed_rhs_rhs1_cpu(
                        typename decltype(derivs_domain)::discrete_element_type(ii, jj))
                        = evaluator.deriv(xN<I1>(), xN<I2>(), ii + shift - 1, jj + shift - 1);
            }
        }
        ddc::Chunk Sderiv_mixed_lhs_lhs1_alloc(
                derivs_domain,
                ddc::KokkosAllocator<double, MemorySpace>());
        ddc::ChunkSpan Sderiv_mixed_lhs_lhs1 = Sderiv_mixed_lhs_lhs1_alloc.span_view();
        ddc::parallel_deepcopy(Sderiv_mixed_lhs_lhs1, Sderiv_mixed_lhs_lhs1_cpu);
        ddc::Chunk Sderiv_mixed_rhs_lhs1_alloc(
                derivs_domain,
                ddc::KokkosAllocator<double, MemorySpace>());
        ddc::ChunkSpan Sderiv_mixed_rhs_lhs1 = Sderiv_mixed_rhs_lhs1_alloc.span_view();
        ddc::parallel_deepcopy(Sderiv_mixed_rhs_lhs1, Sderiv_mixed_rhs_lhs1_cpu);
        ddc::Chunk Sderiv_mixed_lhs_rhs1_alloc(
                derivs_domain,
                ddc::KokkosAllocator<double, MemorySpace>());
        ddc::ChunkSpan Sderiv_mixed_lhs_rhs1 = Sderiv_mixed_lhs_rhs1_alloc.span_view();
        ddc::parallel_deepcopy(Sderiv_mixed_lhs_rhs1, Sderiv_mixed_lhs_rhs1_cpu);
        ddc::Chunk Sderiv_mixed_rhs_rhs1_alloc(
                derivs_domain,
                ddc::KokkosAllocator<double, MemorySpace>());
        ddc::ChunkSpan Sderiv_mixed_rhs_rhs1 = Sderiv_mixed_rhs_rhs1_alloc.span_view();
        ddc::parallel_deepcopy(Sderiv_mixed_rhs_rhs1, Sderiv_mixed_rhs_rhs1_cpu);

        ddc::parallel_for_each(
                exec_space,
                dom_derivs,
                KOKKOS_LAMBDA(typename decltype(dom_derivs)::discrete_element_type const e) {
                    Sderiv_mixed_lhs_lhs(e)
                            = Sderiv_mixed_lhs_lhs1(ddc::select<ddc::Deriv<I1>, ddc::Deriv<I2>>(e));
                    Sderiv_mixed_rhs_lhs(e)
                            = Sderiv_mixed_rhs_lhs1(ddc::select<ddc::Deriv<I1>, ddc::Deriv<I2>>(e));
                    Sderiv_mixed_lhs_rhs(e)
                            = Sderiv_mixed_lhs_rhs1(ddc::select<ddc::Deriv<I1>, ddc::Deriv<I2>>(e));
                    Sderiv_mixed_rhs_rhs(e)
                            = Sderiv_mixed_rhs_rhs1(ddc::select<ddc::Deriv<I1>, ddc::Deriv<I2>>(e));
                });
    }
#endif

    // Instantiate chunk of spline coefs to receive output of spline_builder
    ddc::Chunk coef_alloc(dom_spline, ddc::KokkosAllocator<double, MemorySpace>());
    ddc::ChunkSpan coef = coef_alloc.span_view();

    // Finally compute the spline by filling `coef`
#if defined(BC_HERMITE)
    spline_builder(
            coef,
            vals.span_cview(),
            std::optional(Sderiv1_lhs.span_cview()),
            std::optional(Sderiv1_rhs.span_cview()),
            std::optional(Sderiv2_lhs.span_cview()),
            std::optional(Sderiv2_rhs.span_cview()),
            std::optional(Sderiv_mixed_lhs_lhs.span_cview()),
            std::optional(Sderiv_mixed_rhs_lhs.span_cview()),
            std::optional(Sderiv_mixed_lhs_rhs.span_cview()),
            std::optional(Sderiv_mixed_rhs_rhs.span_cview()));
#else
    spline_builder(coef, vals.span_cview());
#endif
    // Instantiate a SplineEvaluator over interest dimension and batched along other dimensions
#if defined(BC_PERIODIC)
    ddc::PeriodicExtrapolationRule<I1> extrapolation_rule_1;
    ddc::PeriodicExtrapolationRule<I2> extrapolation_rule_2;
#else
    ddc::NullExtrapolationRule extrapolation_rule_1;
    ddc::NullExtrapolationRule extrapolation_rule_2;
#endif
    ddc::SplineEvaluator2D<
            ExecSpace,
            MemorySpace,
            BSplines<I1>,
            BSplines<I2>,
            IDim<I1, I1, I2>,
            IDim<I2, I1, I2>,
#if defined(BC_PERIODIC)
            ddc::PeriodicExtrapolationRule<I1>,
            ddc::PeriodicExtrapolationRule<I1>,
            ddc::PeriodicExtrapolationRule<I2>,
            ddc::PeriodicExtrapolationRule<I2>,
#else
            ddc::NullExtrapolationRule,
            ddc::NullExtrapolationRule,
            ddc::NullExtrapolationRule,
            ddc::NullExtrapolationRule,
#endif
            IDim<X, I1, I2>...>
            spline_evaluator(
                    extrapolation_rule_1,
                    extrapolation_rule_1,
                    extrapolation_rule_2,
                    extrapolation_rule_2);

    // Instantiate chunk of coordinates of dom_interpolation
    ddc::Chunk coords_eval_alloc(dom_vals, ddc::KokkosAllocator<Coord<X...>, MemorySpace>());
    ddc::ChunkSpan coords_eval = coords_eval_alloc.span_view();
    ddc::parallel_for_each(
            exec_space,
            coords_eval.domain(),
            KOKKOS_LAMBDA(Index<IDim<X, I1, I2>...> const e) {
                coords_eval(e) = ddc::coordinate(e);
            });


    // Instantiate chunks to receive outputs of spline_evaluator
    ddc::Chunk spline_eval_alloc(dom_vals, ddc::KokkosAllocator<double, MemorySpace>());
    ddc::ChunkSpan spline_eval = spline_eval_alloc.span_view();
    ddc::Chunk spline_eval_deriv1_alloc(dom_vals, ddc::KokkosAllocator<double, MemorySpace>());
    ddc::ChunkSpan spline_eval_deriv1 = spline_eval_deriv1_alloc.span_view();
    ddc::Chunk spline_eval_deriv2_alloc(dom_vals, ddc::KokkosAllocator<double, MemorySpace>());
    ddc::ChunkSpan spline_eval_deriv2 = spline_eval_deriv2_alloc.span_view();
    ddc::Chunk spline_eval_deriv12_alloc(dom_vals, ddc::KokkosAllocator<double, MemorySpace>());
    ddc::ChunkSpan spline_eval_deriv12 = spline_eval_deriv12_alloc.span_view();

    // Call spline_evaluator on the same mesh we started with
    spline_evaluator(spline_eval, coords_eval.span_cview(), coef.span_cview());
    spline_evaluator
            .template deriv<I1>(spline_eval_deriv1, coords_eval.span_cview(), coef.span_cview());
    spline_evaluator
            .template deriv<I2>(spline_eval_deriv2, coords_eval.span_cview(), coef.span_cview());
    spline_evaluator.template deriv2<
            I1,
            I2>(spline_eval_deriv12, coords_eval.span_cview(), coef.span_cview());

    // Checking errors (we recover the initial values)
    double max_norm_error = ddc::parallel_transform_reduce(
            exec_space,
            spline_eval.domain(),
            0.,
            ddc::reducer::max<double>(),
            KOKKOS_LAMBDA(Index<IDim<X, I1, I2>...> const e) {
                return Kokkos::abs(spline_eval(e) - vals(e));
            });
    double max_norm_error_diff1 = ddc::parallel_transform_reduce(
            exec_space,
            spline_eval_deriv1.domain(),
            0.,
            ddc::reducer::max<double>(),
            KOKKOS_LAMBDA(Index<IDim<X, I1, I2>...> const e) {
                Coord<I1> const x = ddc::coordinate(ddc::select<IDim<I1, I1, I2>>(e));
                Coord<I2> const y = ddc::coordinate(ddc::select<IDim<I2, I1, I2>>(e));
                return Kokkos::abs(spline_eval_deriv1(e) - evaluator.deriv(x, y, 1, 0));
            });
    double max_norm_error_diff2 = ddc::parallel_transform_reduce(
            exec_space,
            spline_eval_deriv2.domain(),
            0.,
            ddc::reducer::max<double>(),
            KOKKOS_LAMBDA(Index<IDim<X, I1, I2>...> const e) {
                Coord<I1> const x = ddc::coordinate(ddc::select<IDim<I1, I1, I2>>(e));
                Coord<I2> const y = ddc::coordinate(ddc::select<IDim<I2, I1, I2>>(e));
                return Kokkos::abs(spline_eval_deriv2(e) - evaluator.deriv(x, y, 0, 1));
            });
    double max_norm_error_diff12 = ddc::parallel_transform_reduce(
            exec_space,
            spline_eval_deriv1.domain(),
            0.,
            ddc::reducer::max<double>(),
            KOKKOS_LAMBDA(Index<IDim<X, I1, I2>...> const e) {
                Coord<I1> const x = ddc::coordinate(ddc::select<IDim<I1, I1, I2>>(e));
                Coord<I2> const y = ddc::coordinate(ddc::select<IDim<I2, I1, I2>>(e));
                return Kokkos::abs(spline_eval_deriv12(e) - evaluator.deriv(x, y, 1, 1));
            });


    double const max_norm = evaluator.max_norm();
    double const max_norm_diff1 = evaluator.max_norm(1, 0);
    double const max_norm_diff2 = evaluator.max_norm(0, 1);
    double const max_norm_diff12 = evaluator.max_norm(1, 1);

    SplineErrorBounds<evaluator_type<IDim<I1, I1, I2>, IDim<I2, I1, I2>>> error_bounds(evaluator);
    EXPECT_LE(
            max_norm_error,
            std::
                    max(error_bounds
                                .error_bound(dx<I1>(ncells), dx<I2>(ncells), s_degree, s_degree),
                        1.0e-14 * max_norm));
    EXPECT_LE(
            max_norm_error_diff1,
            std::
                    max(error_bounds.error_bound_on_deriv_1(
                                dx<I1>(ncells),
                                dx<I2>(ncells),
                                s_degree,
                                s_degree),
                        1e-12 * max_norm_diff1));
    EXPECT_LE(
            max_norm_error_diff2,
            std::
                    max(error_bounds.error_bound_on_deriv_2(
                                dx<I1>(ncells),
                                dx<I2>(ncells),
                                s_degree,
                                s_degree),
                        1e-12 * max_norm_diff2));
    EXPECT_LE(
            max_norm_error_diff12,
            std::
                    max(error_bounds.error_bound_on_deriv_12(
                                dx<I1>(ncells),
                                dx<I2>(ncells),
                                s_degree,
                                s_degree),
                        1e-12 * max_norm_diff12));
}

#if defined(BC_PERIODIC) && defined(BSPLINES_TYPE_UNIFORM)
#define SUFFIX(name) name##Periodic##Uniform
#elif defined(BC_PERIODIC) && defined(BSPLINES_TYPE_NON_UNIFORM)
#define SUFFIX(name) name##Periodic##NonUniform
#elif defined(BC_GREVILLE) && defined(BSPLINES_TYPE_UNIFORM)
#define SUFFIX(name) name##Greville##Uniform
#elif defined(BC_GREVILLE) && defined(BSPLINES_TYPE_NON_UNIFORM)
#define SUFFIX(name) name##Greville##NonUniform
#elif defined(BC_HERMITE) && defined(BSPLINES_TYPE_UNIFORM)
#define SUFFIX(name) name##Hermite##Uniform
#elif defined(BC_HERMITE) && defined(BSPLINES_TYPE_NON_UNIFORM)
#define SUFFIX(name) name##Hermite##NonUniform
#endif

TEST(SUFFIX(Batched2dSplineHost), 2DXY)
{
    Batched2dSplineTest<
            Kokkos::DefaultHostExecutionSpace,
            Kokkos::DefaultHostExecutionSpace::memory_space,
            DimX,
            DimY,
            DimX,
            DimY>();
}

TEST(SUFFIX(Batched2dSplineDevice), 2DXY)
{
    Batched2dSplineTest<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            DimX,
            DimY,
            DimX,
            DimY>();
}

TEST(SUFFIX(Batched2dSplineHost), 3DXY)
{
    Batched2dSplineTest<
            Kokkos::DefaultHostExecutionSpace,
            Kokkos::DefaultHostExecutionSpace::memory_space,
            DimX,
            DimY,
            DimX,
            DimY,
            DimZ>();
}

TEST(SUFFIX(Batched2dSplineHost), 3DXZ)
{
    Batched2dSplineTest<
            Kokkos::DefaultHostExecutionSpace,
            Kokkos::DefaultHostExecutionSpace::memory_space,
            DimX,
            DimZ,
            DimX,
            DimY,
            DimZ>();
}

TEST(SUFFIX(Batched2dSplineHost), 3DYZ)
{
    Batched2dSplineTest<
            Kokkos::DefaultHostExecutionSpace,
            Kokkos::DefaultHostExecutionSpace::memory_space,
            DimY,
            DimZ,
            DimX,
            DimY,
            DimZ>();
}

TEST(SUFFIX(Batched2dSplineDevice), 3DXY)
{
    Batched2dSplineTest<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            DimX,
            DimY,
            DimX,
            DimY,
            DimZ>();
}

TEST(SUFFIX(Batched2dSplineDevice), 3DXZ)
{
    Batched2dSplineTest<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            DimX,
            DimZ,
            DimX,
            DimY,
            DimZ>();
}

TEST(SUFFIX(Batched2dSplineDevice), 3DYZ)
{
    Batched2dSplineTest<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            DimY,
            DimZ,
            DimX,
            DimY,
            DimZ>();
}
