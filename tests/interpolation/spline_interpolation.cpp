// SPDX-License-Identifier: MIT
#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <gtest/gtest.h>

#include "ddc_aliases.hpp"
#include "i_interpolation.hpp"
#include "i_interpolation_builder.hpp"
#include "i_interpolation_evaluator.hpp"
#include "spline_interpolation.hpp"

namespace {

struct X
{
    static bool constexpr PERIODIC = false;
};
struct Y
{
    static bool constexpr PERIODIC = true;
};
struct Z
{
    static bool constexpr PERIODIC = false;
};

int constexpr Degree = 3;

struct BSplinesX : ddc::UniformBSplines<X, Degree>
{
};
struct BSplinesY : ddc::UniformBSplines<Y, Degree>
{
};
struct BSplinesZ : ddc::UniformBSplines<Z, Degree>
{
};

using SplineInterpPointsX = ddc::GrevilleInterpolationPoints<
        BSplinesX,
        ddc::SplineBuilderClosure::GREVILLE,
        ddc::SplineBuilderClosure::GREVILLE>;
using SplineInterpPointsY = ddc::GrevilleInterpolationPoints<
        BSplinesY,
        ddc::SplineBuilderClosure::PERIODIC,
        ddc::SplineBuilderClosure::PERIODIC>;
using SplineInterpPointsZ = ddc::GrevilleInterpolationPoints<
        BSplinesZ,
        ddc::SplineBuilderClosure::GREVILLE,
        ddc::SplineBuilderClosure::GREVILLE>;

struct GridX : SplineInterpPointsX::interpolation_discrete_dimension_type
{
};
struct GridY : SplineInterpPointsY::interpolation_discrete_dimension_type
{
};
struct GridZ : SplineInterpPointsZ::interpolation_discrete_dimension_type
{
};

using ExecSpace = Kokkos::DefaultExecutionSpace;

using Interp1D = detail::SplineInterpolator<
        ExecSpace,
        BSplinesX,
        GridX,
        ddc::detail::TypeSeq<ddc::NullExtrapolationRule, ddc::NullExtrapolationRule>,
        SplineBoundaryClosure::Greville_Greville,
        ddc::SplineSolver::LAPACK>;

using Interp2D = detail::SplineInterpolator2D<
        ExecSpace,
        BSplinesX,
        BSplinesY,
        GridX,
        GridY,
        ddc::detail::TypeSeq<
                ddc::detail::TypeSeq<ddc::NullExtrapolationRule, ddc::NullExtrapolationRule>,
                ddc::detail::
                        TypeSeq<ddc::PeriodicExtrapolationRule<Y>, ddc::PeriodicExtrapolationRule<Y>>>,
        ddc::detail::TypeSeq<
                SplineBoundaryClosure::Greville_Greville,
                SplineBoundaryClosure::Periodic>,
        ddc::SplineSolver::LAPACK>;

using Interp3D = detail::SplineInterpolator3D<
        ExecSpace,
        BSplinesX,
        BSplinesY,
        BSplinesZ,
        GridX,
        GridY,
        GridZ,
        ddc::detail::TypeSeq<
                ddc::detail::TypeSeq<ddc::NullExtrapolationRule, ddc::NullExtrapolationRule>,
                ddc::detail::
                        TypeSeq<ddc::PeriodicExtrapolationRule<Y>, ddc::PeriodicExtrapolationRule<Y>>,
                ddc::detail::TypeSeq<ddc::NullExtrapolationRule, ddc::NullExtrapolationRule>>,
        ddc::detail::TypeSeq<
                SplineBoundaryClosure::Greville_Greville,
                SplineBoundaryClosure::Periodic,
                SplineBoundaryClosure::Greville_Greville>,
        ddc::SplineSolver::LAPACK>;

} // namespace

TEST(SplineInterpolatorAPI, OneD)
{
    static_assert(concepts::Interpolation<Interp1D>);
    static_assert(concepts::Interpolation1D<Interp1D>);
    static_assert(concepts::InterpolationBuilder<Interp1D::BuilderType>);
    static_assert(concepts::InterpolationBuilder1D<Interp1D::BuilderType>);
    static_assert(concepts::InterpolationEvaluator<Interp1D::EvaluatorType>);
    static_assert(Interp1D::rank() == 1);
    static_assert(std::is_same_v<
                   Interp1D::BuilderType,
                   ddc::SplineBuilder<
                           ExecSpace,
                           ExecSpace::memory_space,
                           BSplinesX,
                           GridX,
                           ddc::SplineBuilderClosure::GREVILLE,
                           ddc::SplineBuilderClosure::GREVILLE,
                           ddc::SplineSolver::LAPACK>>);
}

TEST(SplineInterpolatorAPI, TwoD)
{
    static_assert(concepts::Interpolation<Interp2D>);
    static_assert(concepts::InterpolationBuilder<Interp2D::BuilderType>);
    static_assert(concepts::InterpolationEvaluator<Interp2D::EvaluatorType>);
    static_assert(Interp2D::rank() == 2);
}

TEST(SplineInterpolatorAPI, ThreeD)
{
    static_assert(concepts::Interpolation<Interp3D>);
    static_assert(concepts::InterpolationBuilder<Interp3D::BuilderType>);
    static_assert(concepts::InterpolationEvaluator<Interp3D::EvaluatorType>);
    static_assert(Interp3D::rank() == 3);
}

TEST(SplineInterpolatorTraits, RankAndCoeffIdxRangeShape)
{
    static_assert(InterpolationBuilderTraits<Interp1D::BuilderType>::rank() == 1);
    static_assert(InterpolationBuilderTraits<Interp2D::BuilderType>::rank() == 2);
    static_assert(InterpolationBuilderTraits<Interp3D::BuilderType>::rank() == 3);

    static_assert(InterpolationEvaluatorTraits<Interp1D::EvaluatorType>::rank() == 1);
    static_assert(InterpolationEvaluatorTraits<Interp2D::EvaluatorType>::rank() == 2);
    static_assert(InterpolationEvaluatorTraits<Interp3D::EvaluatorType>::rank() == 3);

    static_assert(std::is_same_v<
                   InterpolationBuilderTraits<Interp3D::BuilderType>::coeff_idx_range_type,
                   IdxRange<BSplinesX, BSplinesY, BSplinesZ>>);
}

TEST(SplineInterpolatorResolver, OneD)
{
    using Resolved = SplineInterpolator<
            ExecSpace,
            IdxRange<BSplinesX>,
            IdxRange<GridX>,
            ExtrapolationRule::Null_Null,
            SplineBoundaryClosure::Greville_Greville>;
    static_assert(std::is_same_v<Resolved, Interp1D>);
}

TEST(SplineInterpolatorResolver, TwoD)
{
    using Resolved = SplineInterpolator<
            ExecSpace,
            IdxRange<BSplinesX, BSplinesY>,
            IdxRange<GridX, GridY>,
            ExtrapolationRule::Null_Null,
            ExtrapolationRule::Periodic,
            SplineBoundaryClosure::Greville_Greville,
            SplineBoundaryClosure::Periodic>;
    static_assert(std::is_same_v<Resolved, Interp2D>);
}

TEST(SplineInterpolatorResolver, ThreeD)
{
    using Resolved = SplineInterpolator<
            ExecSpace,
            IdxRange<BSplinesX, BSplinesY, BSplinesZ>,
            IdxRange<GridX, GridY, GridZ>,
            ExtrapolationRule::Null_Null,
            ExtrapolationRule::Periodic,
            ExtrapolationRule::Null_Null,
            SplineBoundaryClosure::Greville_Greville,
            SplineBoundaryClosure::Periodic,
            SplineBoundaryClosure::Greville_Greville>;
    static_assert(std::is_same_v<Resolved, Interp3D>);
}
