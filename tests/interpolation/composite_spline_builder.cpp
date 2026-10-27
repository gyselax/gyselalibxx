// SPDX-License-Identifier: MIT
#include <type_traits>

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <gtest/gtest.h>

#include "composite_spline_builder.hpp"
#include "ddc_aliases.hpp"
#include "i_interpolation_builder.hpp"

namespace {

struct X
{
    static bool constexpr PERIODIC = false;
};
struct Y
{
    static bool constexpr PERIODIC = true;
};
struct Batch
{
};

int constexpr Degree = 3;

struct BSplinesX : ddc::UniformBSplines<X, Degree, X::PERIODIC>
{
};
struct BSplinesY : ddc::UniformBSplines<Y, Degree, Y::PERIODIC>
{
};

auto constexpr ClosureX = ddc::SplineBuilderClosure::GREVILLE;
auto constexpr ClosureY = ddc::SplineBuilderClosure::PERIODIC;

using SplineInterpPointsX = ddc::GrevilleInterpolationPoints<BSplinesX, ClosureX, ClosureX>;
using SplineInterpPointsY = ddc::GrevilleInterpolationPoints<BSplinesY, ClosureY, ClosureY>;

struct GridX : SplineInterpPointsX::interpolation_discrete_dimension_type
{
};
struct GridY : SplineInterpPointsY::interpolation_discrete_dimension_type
{
};
struct GridBatch : UniformGridBase<Batch>
{
};

using IdxRangeXY = IdxRange<GridX, GridY>;
using IdxRangeXBatchY = IdxRange<GridX, GridBatch, GridY>;

using ExecSpace = Kokkos::DefaultExecutionSpace;
using MemSpace = ExecSpace::memory_space;

using BuilderX = ddc::SplineBuilder<
        ExecSpace,
        MemSpace,
        BSplinesX,
        GridX,
        ClosureX,
        ClosureX,
        ddc::SplineSolver::LAPACK>;
using BuilderY = ddc::SplineBuilder<
        ExecSpace,
        MemSpace,
        BSplinesY,
        GridY,
        ClosureY,
        ClosureY,
        ddc::SplineSolver::LAPACK>;

using BuilderXY = CompositeSplineBuilder2D<BuilderX, BuilderY>;
using RefBuilderXY = ddc::SplineBuilder2D<
        ExecSpace,
        MemSpace,
        BSplinesX,
        BSplinesY,
        GridX,
        GridY,
        ClosureX,
        ClosureX,
        ClosureY,
        ClosureY,
        ddc::SplineSolver::LAPACK>;

void init_discrete_spaces()
{
    ddc::init_discrete_space<BSplinesX>(Coord<X>(0.0), Coord<X>(1.0), IdxStep<GridX>(10));
    ddc::init_discrete_space<BSplinesY>(Coord<Y>(0.0), Coord<Y>(2.0), IdxStep<GridY>(12));
    ddc::init_discrete_space<GridX>(SplineInterpPointsX::get_sampling<GridX>());
    ddc::init_discrete_space<GridY>(SplineInterpPointsY::get_sampling<GridY>());
    ddc::init_discrete_space<GridBatch>(
            GridBatch::init(Coord<Batch>(0.0), Coord<Batch>(1.0), IdxStep<GridBatch>(4)));
}

IdxRangeXY get_interpolation_idx_range()
{
    return IdxRangeXY(
            SplineInterpPointsX::get_domain<GridX>(),
            SplineInterpPointsY::get_domain<GridY>());
}

} // namespace

TEST(CompositeSplineBuilder2D, Types)
{
    using RefTraits = InterpolationBuilderTraits<RefBuilderXY>;
    static_assert(std::is_same_v<BuilderXY::exec_space, RefBuilderXY::exec_space>);
    static_assert(std::is_same_v<BuilderXY::memory_space, RefBuilderXY::memory_space>);
    static_assert(std::is_same_v<BuilderXY::data_type, RefTraits::data_type>);
    static_assert(std::is_same_v<
                  BuilderXY::interpolation_idx_range_type,
                  RefTraits::interpolation_idx_range_type>);
    static_assert(std::is_same_v<BuilderXY::coeff_idx_range_type, RefTraits::coeff_idx_range_type>);
    static_assert(std::is_same_v<
                  BuilderXY::batched_basis_idx_range_type<IdxRangeXBatchY>,
                  RefBuilderXY::batched_spline_domain_type<IdxRangeXBatchY>>);
    static_assert(std::is_same_v<
                  BuilderXY::batched_derivs_idx_range_type1<IdxRangeXBatchY>,
                  RefBuilderXY::batched_derivs_domain_type1<IdxRangeXBatchY>>);
    static_assert(std::is_same_v<
                  BuilderXY::batched_derivs_idx_range_type2<IdxRangeXBatchY>,
                  RefBuilderXY::batched_derivs_domain_type2<IdxRangeXBatchY>>);
    static_assert(std::is_same_v<
                  BuilderXY::batched_derivs_idx_range_type<IdxRangeXBatchY>,
                  RefBuilderXY::batched_derivs_domain_type<IdxRangeXBatchY>>);
}

TEST(CompositeSplineBuilder2D, BatchedBasisIdxRange)
{
    init_discrete_spaces();
    IdxRangeXY idx_range_xy = get_interpolation_idx_range();
    IdxRangeXBatchY idx_range_xby(
            ddc::select<GridX>(idx_range_xy),
            IdxRange<GridBatch>(Idx<GridBatch>(0), IdxStep<GridBatch>(4)),
            ddc::select<GridY>(idx_range_xy));

    BuilderX builder_x(idx_range_xy);
    BuilderY builder_y(idx_range_xy);
    BuilderXY builder(builder_x, builder_y);
    RefBuilderXY ref_builder(idx_range_xy);

    EXPECT_EQ(
            builder.batched_basis_idx_range(idx_range_xy),
            ref_builder.batched_spline_domain(idx_range_xy));
    EXPECT_EQ(
            builder.batched_basis_idx_range(idx_range_xby),
            ref_builder.batched_spline_domain(idx_range_xby));
}
